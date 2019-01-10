#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
## @file   test/run.py
## @brief  Script to run unit tests in mini
## @author Sergey Lyskov
from __future__ import print_function

import sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

import os, re, subprocess, time
from os import path
from optparse import OptionParser

import json

#The factor by which to multiply the single test timeout value for a suite
SUITE_FACTOR = 3

UnitTestExecutable = ["protocols.test", "core.test", "basic.test", "ObjexxFCL.test", "numeric.test", "utility.test", "apps.test", "devel.test"]
#UnitTestExecutable = ["apps.test", "numeric.test"]

def execute(message, command_line, return_='status', until_successes=False, terminate_on_failure=True, silent=False, silence_output=False, silence_output_on_errors=False, add_message_and_command_line_to_output=False):
    if not silent: print(message);  print(command_line); sys.stdout.flush();
    while True:

        p = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, errors = p.communicate()
        output = output + errors
        output = output.decode(encoding='utf-8', errors='backslashreplace')
        exit_code = p.returncode

        #exit_code, output = subprocess.getstatusoutput(command_line)

        if (exit_code  and  not silence_output_on_errors) or  not (silent or silence_output): print(output); sys.stdout.flush();

        if exit_code and until_successes: pass  # Thats right - redability COUNT!
        else: break

        print( "Error while executing {}: {}\n".format(message, output) )
        print("Sleeping 60s... then I will retry...")
        sys.stdout.flush();
        time.sleep(60)

    if add_message_and_command_line_to_output: output = message + '\nCommand line: ' + command_line + '\n' + output

    if return_ == 'tuple': return(exit_code, output)

    if exit_code and terminate_on_failure:
        print("\nEncounter error while executing: " + command_line)
        if return_==True: return True
        else:
            print('\nEncounter error while executing: ' + command_line + '\n' + output);
            raise BenchmarkError('\nEncounter error while executing: ' + command_line + '\n' + output)

    if return_ == 'output': return output
    else: return exit_code



# Output info class, simple wrapper for named collection of fields.
class OI:
    def __init__(self, **entries): self.__dict__.update(entries)


class Tester:
    def __init__(self):
        #self.db_path = db_path
        self.systemLog = ""  # System log - we store all information here
        self.results = {}
        self.jobs = []  # list of spawned process pid's
        self.platform = None
        self.testpath = None


    # print and log given mesage, return message
    def log(self, s):
        self.systemLog += s
        print(s)
        return s


    def mfork(self):
        ''' Check if number of child process is below Options.jobs. And if it is - fork the new pocees and return its pid.
        '''
        while len(self.jobs) >= Options.jobs :
            for p in self.jobs[:] :
                r = os.waitpid(p, os.WNOHANG)
                if r == (p, 0):  # process have ended without error
                    self.jobs.remove(p)
                elif r[0] == p :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                    for p in self.jobs: os.waitpid(p, 0)
                    print('Some of the unit test suite terminate abnormally!')
                    sys.exit(1)

            if len(self.jobs) >= Options.jobs: time.sleep(.5)
        pid = os.fork()
        if pid: self.jobs.append(pid) # We are parent!
        return pid


    # Try to identity plaform by using scons compiliation feature.
    def getPlatformID(self):
        if Options.CMake:
            self.testpath = Options.testpath
            self.log( "Skipping platform identification. Using explicit directory instead: " + self.testpath )
            return
        self.log( "Identifying platform...\n")
        cmd_str = "./scons.py unit_test_platform_only log=platform"
        if Options.extras:
            cmd_str += " extras=%s" % (Options.extras)
        if Options.mode:
            cmd_str += " mode=%s" % (Options.mode)
        if Options.compiler:
            cmd_str += " cxx=%s" % (Options.compiler)
        if Options.cxx_ver:
            cmd_str += " cxx_ver=%s" % (Options.cxx_ver)

        pl = subprocess.check_output(cmd_str, shell=True).decode('utf-8', errors="replace")
        lines = pl.split('\n')
        for s in lines:
            if  len( s.split() ) > 1 and s.split()[0] == 'Platform:':
                platform = s.split()[1]
                self.log( "Platform found: " + platform )
                self.platform = platform
                self.testpath = "build/test/" + platform
                return
        sys.exit("run.py is about to crash because it could not use SCons to detect your platform.  The most likely reason for this is that you are running it from the wrong directory - it must be run from the rosetta_source directory, not the rosetta_source/test directory, even though it lives in the latter.\nScons output for command line {0} is:\n{1}\n".format(cmd_str, pl))
        return "PlatformWasNotFound!!!"  # <-- That should not really happend.


    # extract information regarding how good unit tests run.
    def extractInfo(self, output):
        # default init, in case we can't extract information
        oi = OI(testCount=0, testFailed=0, failedTestsList=[])

        # extracting number of test
        s = re.search(r'Running (\d\d*) test', output)
        if s:
            g = s.groups()
            oi.testCount = int( g[0] )

        # extracting number of test failed
        s = re.search(r'Failed (\d\d*) of \d\d* test', output)
        if s:
            g = s.groups()
            oi.testFailed = int( g[0] )

            # extracting names of the failed tests
            s = re.findall(r'CXXTEST_ERROR: (\w*) Failed!', output)
            for t in s: oi.failedTestsList.append(t)

            #s = re.search(r'CXXTEST_ERROR: (\w*) Failed!', output)
            #if s:
            #    g = s.groups()
            #    for t in g: oi.failedTestsList.append(t)

        else: # all test pussed?
            s = re.search(r'All tests passed!', output)
            if not s: # Something wrong then, count all tests as failed.
                oi.testFailed = oi.testCount
                oi.failedTestsList = ['All']

        #print oi.__dict__
        return oi


    def runOneLibUnitTests(self, lib, yjson_file, log_file):
        if Options.one  and  ( (lib, Options.one.split(':')[0]) not in self.all_test_suites_by_lib[lib] ): return

        #self.unitTestLog += self.log("-------- %s --------\n" % E)
        path = "cd " + self.testpath + " && "
        exe = self.buildCommandLine(lib, Options.one)
        print("Command line:: %s%s" % (path, exe))

        if os.path.isfile(yjson_file): os.remove(yjson_file)

        output = [ "Running %s unit tests...\n" % lib ]

        if( Options.timeout > 0 ):
            timelimit = sys.executable + ' ' +os.path.abspath('test/timelimit.py') + ' ' + str(Options.timeout) + ' '
        else:
            timelimit = ""
        #print "timelimit:", timelimit

        command_line = path + ' ' + timelimit + exe + " 1>&2"
        #print 'Executing:', command_line

        f = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr
        if Options.valgrind:
            self.parseOneLibValgrind(f, output, yjson_file)
        else:
            self.parseOneLib(f, output)

        f.close()

        # saving log to a file...
        with open(log_file, 'w') as f:
            f.write(''.join(output))

    def parseOneLib(self, f, output):
        '''Parse output for OneLib run for regular unit tests.'''
        for line in f:
            l = line.decode('utf-8', errors="replace")
            print(l, end='')
            output.append(l)
            sys.stdout.flush()

    def parseOneLibValgrind(self, f, output, yjson_file):
        valgrind_errors = 0
        for line in f:
            if line.startswith("=="):
                print(line.strip())
                if line.find("ERROR SUMMARY") != -1:
                    valgrind_errors += int( line.split()[3] )
                output.append(line)
                sys.stdout.flush()

        #The yjson file might not be the same name - e.g. if we're running a single suite
        if valgrind_errors > 0 and os.path.isfile(yjson_file):
            data = json.loads(open(yjson_file).read())
            failures = []
            for test in data["ALL_TESTS"]:
                if test.startswith( Options.one ):
                    failures.append( test )
            data["FAILED_TESTS"] = failures
            yjson = open( yjson_file, "w" )
            yjson.write( json.dumps(data) )
            yjson.close()

    def runOneSuite(self, lib, suite):
        path = "cd " + self.testpath + " && "
        exe = self.buildCommandLine(lib, suite)

        log_file = self.testpath + '/' + lib + '.' + suite + '.log'
        yjson_file = self.testpath + '/' + lib + '.'+ suite + '.json'

        if os.path.isfile(yjson_file): os.remove(yjson_file)
        if os.path.isfile(log_file): os.remove(log_file)

        #output = [ ]
        output = [ "Running %s:%s unit tests...\n" % (lib, suite) ]
        #print output

        if Options.timeout > 0:
            timelimit = sys.executable + ' ' + os.path.abspath('test/timelimit.py') + ' ' + str(Options.timeout*SUITE_FACTOR) + ' '
        else:
            timelimit = ""
        #print "timelimit:", timelimit

        command_line = path + ' ' + timelimit + exe + " 1>&2"
        #print 'Executing:', command_line
        f = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr
        if Options.valgrind:
            self.parseOneSuiteValgrind(f, output, lib, suite, yjson_file)
        else:
            self.parseOneSuite(f, output)
        f.close()

        # saving log to a file...
        with open(log_file, 'w') as f:
            f.write(''.join(output))


    def parseOneSuite(self, f, output):
        for line in f:
            l = line.decode('utf-8', errors="replace")
            if l != 'All tests passed!\n':
                print(l, end='')
                output.append(l)
                sys.stdout.flush()

    def parseOneSuiteValgrind(self, f, output, lib, suite, yjson_file):
        valgrind_errors = 0
        for line in f:
            if line != 'All tests passed!\n':
                if line.find("ERROR SUMMARY") != -1:
                    valgrind_errors += int( line.split()[3] )
                    if valgrind_errors != 0:
                        print(lib, suite, line.strip())
                output.append(line)
                sys.stdout.flush()

        # Running suites, we really can't determine which of the sub-tests caused the Valgrind error - so mark them all failed.
        if valgrind_errors > 0:
            data = json.loads(open(yjson_file).read())
            failures = []
            for test in data["ALL_TESTS"]:
                if test.startswith( suite ):
                    failures.append( test )
            data["FAILED_TESTS"] = failures
            yjson = open( yjson_file, "w" )
            yjson.write( json.dumps(data) )
            yjson.close()

    # Run unit test.
    def runUnitTests(self):
        self.getPlatformID()
        #self.log( "Run unit tests...\n")
        self.unitTestLog = "================================ UnitTest Results ================================\n"

        logs_yjsons = {}

        # Getting list of test suites
        self.all_test_suites = []
        self.all_test_suites_by_lib = {}

        self.all_tests = []

        for lib in UnitTestExecutable:
            tests = set()

            exit_code, command_output = execute('Getting list of test for {lib}...'.format(**locals()), self.testpath + '/' + lib + ' _ListAllTests_', return_='tuple', silent=True)

            if exit_code:
                # The output doesn't contain a test list, only an error message
                print('Could not acquire list of available tests from library {lib}, request terminated with error:\n{command_output}\nTerminating...'.format(**locals()))
                sys.exit(1)
            else:
                for test in command_output.split():
                    tests.add( (lib, test.split(':')[0]) )
                    self.all_tests.append(test)

            self.all_test_suites.extend( tests )
            self.all_test_suites_by_lib[lib] = tests


        self.all_tests.sort( key = lambda x: x.lower() )
        #print 'Suites:', self.all_test_suites
        #print 'Tests:', self.all_tests

        if Options.one:  # or Options.jobs < 5:
            if Options.one not in [s for (l,s) in self.all_test_suites] + self.all_tests:
                print('\nTest suite %s not found!' % Options.one)
                print( "Available test suites are: {tests}".format(tests = ' '.join(self.all_tests) ) )
                sys.exit(1)

            for lib in UnitTestExecutable:
                log_file = self.testpath + '/' + lib + '.log'
                yjson_file = self.testpath + '/' + lib + '.json'

                logs_yjsons[lib] = (log_file, yjson_file)

                if Options.jobs > 1:
                    pid = self.mfork()
                    if not pid:  # we are child process
                        self.runOneLibUnitTests(lib, yjson_file, log_file)
                        sys.exit(0)

                else:
                    self.runOneLibUnitTests(lib, yjson_file, log_file)

            for p in self.jobs: os.waitpid(p, 0)  # waiting for all child process to termintate...

            # extracting and aggregation results, but only if running all tests
            if not Options.one:
                all_results_yjson_file = '.unit_test_results.json'
                with open(all_results_yjson_file, 'w') as uf:
                    for lib in logs_yjsons:
                        (log_file, yjson_file) = logs_yjsons[lib]
                        with open(log_file, 'r') as f:
                            info = self.extractInfo( f.read() )
                        info.name = lib
                        self.results[lib] = info

                        if not os.path.isfile(yjson_file):
                            print("Unable to read json file with test results %s - unit test run aborted!" % yjson_file)
                            os.remove(all_results_yjson_file)
                            sys.exit(1)

                        # generating summary yjson file for all tests
                        uf.write(lib + ':\n')
                        with open(yjson_file, 'w') as f:
                            for line in f: uf.write("    "+line)
                        uf.write('\n')

            #self.log( "Run unit tests... Done.\n")

        else: # running Unit test on multiple CPU's, new style, fully parallel
            for lib, suite in self.all_test_suites:
                pid = self.mfork()
                if not pid:  # we are child process
                    if Options.debug: print('DEBUG MODE ENABLED: Skipping tests run: ' + lib + suite)
                    else: self.runOneSuite(lib, suite)
                    sys.exit(0)

            for p in self.jobs: os.waitpid(p, 0)  # waiting for all child process to termintate...

            # Now all tests should be finished, all we have to do is to create a log file and aggegated JSON file to emulate single CPU out run
            all_yjson = {}
            all_json = dict(tests={}, summary=dict(total=0, failed=0), config=dict(test_path=self.testpath), runtime={})
            _failed_ = 'failed'
            _finished_ = 'passed'

            for lib in UnitTestExecutable:
                log_file = self.testpath + '/' + lib + '.log'
                yjson_file = self.testpath + '/' + lib + '.json'

                logs_yjsons[lib] = (log_file, yjson_file)

                yjson_data = {}
                with open(log_file, 'w') as log_file_h:
                    for l, suite in self.all_test_suites_by_lib[lib]:
                        suite_log = open(self.testpath + '/' + lib + '.' + suite + '.log').read()
                        log_file_h.write( suite_log )

                        print('trying: ', self.testpath + '/' + lib + '.' + suite + '.json')
                        with open(self.testpath + '/' + lib + '.' + suite + '.json') as f:
                            data = json.loads( f.read() )

                        for k in data:
                            if k in yjson_data:
                                if type(data[k]) == list:
                                    yjson_data[k] = list( set(yjson_data[k] + data[k]) )
                                else: yjson_data[k].update(data[k])

                            else: yjson_data[k] = data[k]

                        def json_key_from_test(test): return lib[:-len('.test')] + ':' + test.replace(':', ':')  # we might want different separator then ':' in the future

                        for test in data['ALL_TESTS']:
                            json_key = json_key_from_test(test)
                            if json_key not in all_json['tests']: all_json['tests'][json_key] = dict(state=_finished_, log='')

                        for test in data['FAILED_TESTS']:
                            json_key = json_key_from_test(test)
                            all_json['tests'][json_key] = dict(state=_failed_, log=suite_log)

                if 'runtime' in yjson_data: all_json['runtime'].update( yjson_data['runtime'] )

                with open(yjson_file, 'w') as yjson_file_h:
                    yjson_file_h.write( json.dumps(yjson_data) )

                #if 'ALL_TESTS' in yjson_data:
                if yjson_data:
                    self.results[lib] = OI(testCount=len(yjson_data['ALL_TESTS']), testFailed=len(yjson_data['FAILED_TESTS']), failedTestsList=yjson_data['FAILED_TESTS'], name=lib)
                    all_yjson[lib] = yjson_data

                logs_yjsons[lib] = (log_file, yjson_file)

            #print 'All_yjson:', all_yjson
            # with open('.unit_test_results.json', 'w') as f:
            #     f.write( json.dumps(all_yjson) )

            for k in list( all_json['runtime'].keys() ): all_json['runtime'][ k.replace('./', '').replace('.test', '') ] = all_json['runtime'].pop(k)

            for t in self.results:
                all_json['summary']['total']   += self.results[t].testCount
                all_json['summary']['failed'] += self.results[t].testFailed

            all_json['summary']['failed_tests'] = [ t for t in all_json['tests'] if all_json['tests'][t]['state'] == _failed_ ]

            with open('.unit_test_results.json', 'w') as f:
                json.dump(all_json, f, sort_keys=True, indent=2)


        '''
        error = False
        for p in self.jobs:
            r = os.waitpid(p, 0)  # waiting for all child process to termintate...
            if r[0] == p :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                error = True

        if error:
            print 'Some of the build scripts return an error, PyRosetta build failed!'
            sys.exit(1)
        '''


    def printSummary(self):
        total = 0
        failed = 0
        failedTestsList = []

        for t in self.results:
            total += self.results[t].testCount
            failed += self.results[t].testFailed
            failedTestsList.extend( [ self.results[t].name + ": " + r  for r in self.results[t].failedTestsList] )
            #print "--> ", self.results[t].failedTestsList

        print("-------- Unit test summary --------")
        print("Total number of tests:", total)
        print("  number tests passed:", total-failed)
        print("  number tests failed:", failed)
        if failedTestsList:
            print("  failed tests:")
            for t in failedTestsList:
                print("   ", t)
        if total == 0:
            print("Success rate: 0%")
        else:
            print("Success rate: %s%%" % ((total-failed)*100/total))
        print("---------- End of Unit test summary")


    def buildCommandLine(self, lib, test):
        if Options.valgrind:
            preamble = [ Options.valgrind_path ]
            if Options.trackorigins:
                preamble.append( "--track-origins=yes" )
            if Options.leakcheck:
                preamble.append( "--leak-check=full" )
        else:
            preamble = []

        return " ".join( preamble + ["./" + lib, test] +  self.genericTestFlags() )

    def genericTestFlags(self):
        """Generate generic command line flags (database/tracer/etc...) for tests."""

        flags = ["--database", Options.database]

        if Options.mute:
            flags.extend(["-mute"] +  Options.mute)

        if Options.unmute:
            flags.extend(["-unmute"] +  Options.unmute)

        if Options.levels:
            flags.extend(["-out:levels"] + Options.levels)

        if Options.level:
            flags.extend(["-out:level"] + Options.level)

        #No personal flag configs
        flags.append("-no_fconfig")

        return flags

def parse_valgrind_options(options, option_parser):
    if( options.timeout == option_parser.get_default_values().timeout ):
        print("Default value for timeout used - turning off timeout for valgrind.") # Valgrind runs take a long time.
        options.timeout = 0
    if options.valgrind_path is None:
        import distutils.spawn
        options.valgrind_path = distutils.spawn.find_executable('valgrind')
        if options.valgrind_path is None:
            print("Unable to find valgrind - install or specify the path with the --valgrind_path option.")
            sys.exit(1)
    options.valgrind_path = path.abspath( options.valgrind_path )
    if not os.path.exists( options.valgrind_path ):
        print("Cannot find Valgrind at", options.valgrind_path, "install or use the --valgrind_path option.")
        sys.exit(1)

def main(args):
    ''' Script to run Unit test in Rosetta3.  For debugging, you must use the --unmute option to unmute tracers you are interested in.
    '''
    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)

    parser.add_option("-d", '--database',
      default="", #processed below
      action="store",
      help="Path to Rosetta database. (default: $ROSETTA3_DB, ../database, ~/database, ../rosetta_database, ~/rosetta_database, )",
    )

    parser.add_option("--extras",
      action="store",
      help="Extras option passed to scons to help identify platform"
    )

    parser.add_option("--mode",
      action="store",
      help="mode option passed to scons to help identify platform"
    )


    parser.add_option("-1", "--one",
      default='', action="store",
      help="Run just one unit test or one test suite. To run one unit test-suite just specify it name, ie: '--one PDB_IO'. \
            To run one unit test specify name of test suite, colon, then name of the test it self, like this: '--one PDB_IO:test_pdb_io'.  \
            It is HIGHLY RECOMMENDED to run one unit test suite instead of one unit test! (when running one test all suites have to be initialized, while when running one suite this is not the case).",
    )

    parser.add_option("-j", "--jobs",
      default=1,
      type="int",
      help="Number of processors to use when running tests (default: 1)",
    )

    parser.add_option('--mute',
      default=[],
      action="append",
      help="Mute specified tracer channels. (Note - ALL are MUTED by default!)",
    )

    parser.add_option('--unmute',
      default=[],
      action="append",
      help="UnMute specified tracer channels. (You will need this for debugging as ALL CHANNELS ARE MUTED BY DEFAULT)",
    )

    parser.add_option('--levels',
      default=[],
      action="append",
      help="Tracer out:levels configuration.",
    )

    parser.add_option('--level',
      action="append",
      help="Tracer out:level configuration (For example, 400/500)",
    )

    parser.add_option("-c", '--compiler',
      default=None,
      action="store",
      help="Name of the compiler used.",
    )

    parser.add_option('--cxx_ver',
      default=None,
      action="store",
      help="Version of the compiler used.",
    )

    parser.add_option("-C", '--CMake',
      default=False,
      action="store_true",
      help="Was CMake used to build the unit tests? (i.e. obey the --testpath option).",
    )

    parser.add_option("-T", '--testpath',
      default="cmake/build_debug/",
      action="store",
      help="The relative directory where the unit tests were built into. (default: 'cmake/build_debug/')",
    )

    parser.add_option("--debug",
      action="store_true", default=False,
      help="Enable debug mode: skip test running, only generate summary from previously created JSON and log files.",
    )

    parser.add_option("--timeout",
      action="store",type="int",default=10,metavar="MINUTES",
      help="Automatically cancel test runs if they are taking longer than the given number of minutes. (Default 10) A value of zero turns off the timeout."
    )

    parser.add_option("--valgrind",
      default=False, action="store_true",
      help="Enable valgrind checking mode, instead of regular unit testing.",
    )
    parser.add_option("--valgrind_path",
      default=None,
      help="The path to the valgrind executable. (default: find in path)",
    )
    parser.add_option("--trackorigins",
      action="store_true",
      help="Tell Valgrind to output information about where uninitialized variables came from (default is no tracking, which is faster)",
    )
    parser.add_option("--leakcheck",
      action="store_true",
      help="Tell Valgrind to output details about memory leaks in addition to finding uninitialized variables. (default is no info, which is faster)",
    )

    (options, args) = parser.parse_args(args=args[1:])

    unknown_args = [a.strip() for a in args if a.strip()]
    if unknown_args:
        raise ValueError("Unknown options: %s", unknown_args)

    if options.database == parser.get_default_values().database:
        if path.isdir( path.join( "..", "database") ):
            options.database = path.abspath(path.join( "..", "database"))

        elif os.environ.get('ROSETTA3_DB') is not None and \
                path.isdir(os.environ.get('ROSETTA3_DB')):
            options.database = os.environ.get('ROSETTA3_DB')
        elif path.isdir( path.join( "..", "database") ):
            options.database = path.abspath(path.join( "..", "database"))
        elif path.isdir( path.join( path.expanduser("~"), "database") ):
            options.database = path.join( path.expanduser("~"), "database")
        elif path.isdir( path.join( "..", "rosetta_database") ):
            options.database = path.abspath(path.join( "..", "rosetta_database"))
        else:
            raise ValueError("Can't find database at %s; please set $ROSETTA3_DB or use -d" % options.database)

    else: options.database = path.abspath( options.database )


    if not options.valgrind and ( options.trackorigins or options.leakcheck or options.valgrind_path is not None ):
        print("Valgrind specific options set, but Valgrind mode not enabled. Enabling it for you.")
        options.valgrind = True
    if options.valgrind:
        parse_valgrind_options(options,parser)

    global Options;  Options = options


    #db_path = None
    #if len(args) <= 1: print "Warning: no database path was specified!"
    #else: db_path = " ".join(args[1:])

    T = Tester()
    T.runUnitTests()
    if not Options.one: T.printSummary()

    print("Done!")


if __name__ == "__main__": main(sys.argv)
