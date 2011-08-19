#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#
# @author Sergey Lyskov
#

import os, sys, os.path, time, commands, subprocess, datetime
from optparse import OptionParser

# Create global 'Platform' that will hold info of current system
if sys.platform.startswith("linux"): Platform = "linux" # can be linux1, linux2, etc
elif sys.platform == "darwin" : Platform = "macos"
else: Platform = "_unknown_"
#PlatformBits = platform.architecture()[0][:2]


def getSubprocessMemoryInfo(parentPID):
    size = 0
    for line in commands.getoutput('ps --ppid %s -o vsize,pid' % parentPID).split('\n')[1:]:
        sz, pid = map(int, line.split())
        size += sz + getSubprocessMemoryInfo(pid)

    return size




def run(test, options):
    print 'Running test %s...' % test

    workdir = os.path.abspath( os.path.join("tests", test) )
    print ' Test working dir is: %s' % workdir

    # Running tests
    platform = Platform
    minidir = os.path.abspath( './../../' );         print '    Mini home dir is: %s' % minidir
    bin = os.path.join(minidir, "bin")
    compiler = 'gcc'
    mode = 'release'
    binext = platform+compiler+mode

    database = os.path.abspath( options.database );  print 'Mini database dir is: %s' % database

    templates = dict(minidir=minidir, database=database, workdir=workdir, platform=platform, bin=bin, compiler=compiler, mode=mode, binext=binext)

    fname = os.path.join(workdir, 'command')
    cmd = file(fname).read().strip()
    cmd = cmd % templates  # variable substitution using Python printf style
    # Writing result to .sh file for future refference.
    f = file(fname+'.sh', 'w');  f.write(cmd);  f.close()

    output_dir = os.path.join(workdir, 'output')
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Now we actualy run test...
    proc = subprocess.Popen(["bash", fname+'.sh'], preexec_fn=os.setpgrp)

    start = time.time() # refined start time
    #if self.timeout == 0:
    # retcode = proc.wait() # does this block all threads?
    timeout = options.timeout
    if timeout == 0 : timeout = 99999999

    memory = [(0, 0)]
    time.sleep(1)

    while time.time() - start <= timeout:
        retcode = proc.poll()
        if retcode is not None: break
        memory.append( (int(time.time() - start), getSubprocessMemoryInfo(os.getpid())/1000. ) )
        #print 'Time: %s, Memory: %s' % memory[-1]
        time.sleep(1)

    # writing results to a file and generating memory graph
    memory_data_fn = output_dir + '/memory.data'
    f = file(memory_data_fn, 'w')
    f.write('%-10s %-10s\n' % ('time', 'memory'));  map(lambda x: f.write('%-10s %-10s\n' % x), memory);  f.close()

    print commands.getoutput("""echo 'set terminal png small;set xlabel \"time sec\";set ylabel \"memory MB\";plot "%s" u 1:2 notitle' | gnuplot > %s/memory.png""" % (memory_data_fn, output_dir))

    max_memory_allocated = max( map(lambda x: x[1], memory) )
    execution_time = memory[-1][0]

    yaml_data = { 'execution_time_s' : execution_time, 'max_memory_allocated_MB': max_memory_allocated }
    f = file(output_dir + '/.results.yaml', 'w');  f.write( str(yaml_data) );  f.close()


    if retcode is None:
        print "*** Test %s exceeded the timeout and will be killed!" % test
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    if retcode != 0 and retcode is not None:
        error_string = "*** Test %s did not run!  Check your --mode flag and paths. [%s]\n" % (test, datetime.datetime.now())
        print error_string,

        # Writing error_string to a file, so integration test should fail for sure
        file(os.path.join(workdir, ".test_did_not_run.log"), 'w').write(error_string)






def main(argv):
    '''A simple system for running protocol profile end-to-end tests in Mini.
    '''

    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)

    parser.add_option("-d", "--database",
      default=os.path.join( os.path.expanduser("~"), "minirosetta_database"),
      help="Directory where Mini database is found (default: ~/minirosetta_database)",
    )

    parser.add_option("-t", "--timeout",
      default=0,
      type="int",
      help="Maximum runtime for each test, in minutes (default: no limit)",
      metavar="MINUTES",
    )


    (options, args) = parser.parse_args(argv)

    if not os.path.isdir(options.database):
        print "Can't find database at %s; please use -d" % options.database
        return 1
    # Normalize path before we change directories!
    options.database = os.path.abspath(options.database)

    # Make sure the current directory is the script directory:
    # Using argv[] here causes problems when people try to run the script as "python integration.py ..."
    #os.chdir( path.dirname(sys.argv[0]) ) # argv[0] is the script name
    if not os.path.isdir("tests"):
        print "You must run this script from mini/test/profile/"
        return 2

    # Each test consists of a directory with a "command" file in it.
    if len(args) > 0:
        tests = args
    else:
        tests = [ d for d in os.listdir("tests") if not d.startswith(".") and os.path.isdir(os.path.join("tests", d)) ]

    # Now actually running the tests...
    for test in tests:
        run(test, options)





if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
