import os, stat, shutil, sys

# Load projects.settings file and put resultant data into PROJECT_SETTINGS
# This script (build_util.py) is in rosetta_source/cmake/, and want rosetta_source/projects.settings
PROJECTS_SETTINGS = {}
execfile(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../projects.settings") ,PROJECTS_SETTINGS)
KNOWN_PROJECTS = PROJECTS_SETTINGS["projects"]["src"]
KNOWN_TESTS =PROJECTS_SETTINGS["projects"]["test"]

def list_project_files(path_to_mini, project_name):
	path_to_project = path_to_mini + 'src/'
	settings_file_name = path_to_project + project_name + '.src.settings'
	# brutal hack for pilot_apps.src.settings.all
	if project_name == 'pilot_apps':
		settings_file_name += '.all'

	settings_dict = {}
	execfile(settings_file_name, settings_dict)

	project_files = []
	for dir, srcfiles in settings_dict['sources'].iteritems():
		if len(dir) != 0 and dir[-1] != '/':
			dir += '/'

		hash_index = project_name.find('#')
		if hash_index == -1:
			project_name_nohash = project_name
		else:
			project_name_nohash = project_name[:hash_index]
		#full_path = path_to_project + '/' + project_name_nohash + '/' + dir
		full_path = path_to_project + '/' + dir

		# another brutal hack for pilot_apps.src.settings.all
		if project_name == 'pilot_apps' or project_name == 'apps':
			full_path = path_to_project + '/apps/' + '/' + dir

		old_srcfiles = srcfiles
		srcfiles = []
		for srcfile in old_srcfiles:
			if len(srcfile) > 3 and srcfile[-3:] == ".cu" :
				continue # APL NOTE: cannot currently compile cuda files with cmake
			if os.path.exists(full_path + srcfile + '.cc'):
				srcfiles.append(srcfile + '.cc')
			elif os.path.exists(full_path + srcfile + '.c'):
				srcfiles.append(srcfile + '.c')
			else:
				raise RuntimeError('Nonexistant source file: ' + full_path + srcfile)

		if not os.path.exists(full_path):
			# Git ignores empty directories, so if there aren't any headers or source files the listdir() below may choke.
			# A missing directory where something was expected should error out above.
			continue

		hdrfiles = os.listdir(full_path)
		hdrfiles = [hdrfile for hdrfile in hdrfiles if hdrfile.endswith('.hh') or hdrfile.endswith('.h')]

		dirfiles = sorted(hdrfiles + srcfiles)
		project_files.append((dir, dirfiles))

	#return (path_to_project + project_name_nohash + '/', project_files)
	return (path_to_project + '/', project_files)

def list_test_files(path_to_mini, test_name):
	path_to_test = path_to_mini + 'test/'
	settings_file_name = path_to_test + test_name + '.test.settings'
	settings_dict = {}
	execfile(settings_file_name, settings_dict)

	test_files = []
	for dir, srcfiles in settings_dict['sources'].iteritems():
		if len(dir) != 0 and dir[-1] != '/':
			dir += '/'
		dir = test_name + '/' + dir

		dirfiles = []
		for file in srcfiles:
			dirfiles.append((file + '.cxxtest.hh', file + '.cxxtest.cc', False))

		hdrfiles = os.listdir(path_to_test + '/' + dir)
		hdrfiles = [hdrfile for hdrfile in hdrfiles if hdrfile.endswith('.hh') and not hdrfile.endswith('.cxxtest.hh') ]
		for file in hdrfiles:
			dirfiles.append((file, None, None))

		test_files.append((dir, dirfiles))

	dirfiles = []
	dirfiles.append((test_name + '.cxxtest.hh', test_name + '.cxxtest.cc', True))

	hdrfiles = os.listdir(path_to_test + '/' + test_name)
	hdrfiles = [hdrfile for hdrfile in hdrfiles if hdrfile.endswith('.hh') and not hdrfile.endswith('.cxxtest.hh') ]
	for file in hdrfiles:
		dirfiles.append((file, None, None))

	test_files.append((test_name + '/', dirfiles))

	test_inputs_dict = {}
	for input in settings_dict['testinputfiles']:
		filename = os.path.basename(input)
		dir = os.path.dirname(input)
		if len(dir) != 0 and dir[-1] != '/':
			dir += '/'
		dir = test_name + '/' + dir

		if test_inputs_dict.has_key(dir):
			test_inputs_dict[dir].append(filename)
		else:
			test_inputs_dict[dir] = [filename]

	test_inputs = []
	for dir, inputfiles in test_inputs_dict.iteritems():
		test_inputs.append((dir, inputfiles))

	return (path_to_test, test_files, test_inputs)

def update_libraries_list(projects):
	exclude = ['apps', 'pilot_apps']
	out = open('build/libraries.cmake', 'w')
	out.write('SET( LIBRARIES ' + reduce(lambda x,y: x + ' ' + y, filter(lambda x: x not in exclude, projects)) + ')')

def update_test_list(projects, test_path):
	out = open('build/test_libraries.cmake', 'w')
	try:
		out.write('SET( TEST_LIBRARIES ' + ' '.join( projects) + ')\n')
		out.write('SET( SRCDIR ../' + test_path + ' )\n' )
	finally:
		out.close()

def test_main(path_to_mini, argv, project_test_callback = None):

	tests =  [ "all" ]

	if tests == [ "all" ]:
		tests = KNOWN_TESTS

	for test in tests:

		if test not in KNOWN_TESTS:
			print 'unknown test: ' + test
			sys.exit(-1)

		test_path, test_files, test_inputs = list_test_files(path_to_mini, test)
		project_test_callback(test, path_to_mini, "build/", test_files, test_inputs)

	update_test_list(tests, path_to_mini + "test/")

def project_main(path_to_mini, argv, project_callback):
	if len(argv) < 2:
		print 'usage: %s [project]...' % argv[0]
		print '  known projects:'
		for p in KNOWN_PROJECTS + ['all']:
			print '    %s' % p
		sys.exit(-1)
	update_libraries_list(KNOWN_PROJECTS)
	projects = argv[1:]
	if projects == ['all']:
		projects = KNOWN_PROJECTS

	for project in projects:
		if project not in KNOWN_PROJECTS:
			print 'unknown project: ' + project
			sys.exit(-1)

		project_path, project_files = list_project_files(path_to_mini, project)
		project_callback(project, project_path, project_files)


