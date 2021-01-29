import os
import sys

from helpers import read_prm



###########################################################################
## EXECUTE CODE AND check if stdout exists
###########################################################################
def execute(exec_path, prm_path, projectname, log = True, ntail = 0 , mpi_procs = 1 ,**kwargs) :

   cmd = ('${MPIRUNCOMMAND} %d ' % mpi_procs) #should be set in environment (export MPIRUNCOMMAND='mpirun -np' or 'srun -p')

   cmd=cmd+' '+exec_path.strip()+' '+prm_path.strip()
   if log :
      #log_path = os.path.splitext(os.path.basename(prm_path))[0] + '.log'
      log_path = 'log_'+projectname+'.txt'
      plines = open(prm_path, 'r').readlines()
      logf = open(log_path, 'w')
      for pline in plines :
         logf.write(pline)
      logf.close()
   print('execute: '+cmd)
   os.system(cmd+' 2>std.err 1>std.out')
   stdout=open('std.out','r').readlines()
   stderr=open('std.err','r').readlines()
   if log :
      logf = open(log_path, 'a')
      logf.write('EXECUTE COMMAND:\n')
      logf.write(cmd.strip() + '\n')
      logf.write('EXECUTE: STDOUT FILE:\n')
      for line in stdout :
         logf.write(line)
      logf.write('EXECUTE: STDERR FILE:\n')
      for line in stderr :
         logf.write(line)
      logf.close()
   finishedline = [ line for line in stdout[-5:] if ('finished' in line.lower())]
   if (len(finishedline) == 0) :
      sys.exit(99)  ## be sure that execution marked "crashed"
   if (ntail > 0) :
      for line in stdout[-ntail:] :
        print( line.split('\n')[0] )

###########################################################################
## BUILD CODE AND EXECUTE GIVEN TESTS
###########################################################################
def test_fluxo(jobname='test',  case=0, stage=0, base_opts={}, build_opts={},force_opts={}, ntail=0, mpi_procs=1, err=None, run_opts=None, **kwargs ):
   # build code with buildopts, 
   cwd = os.getcwd()  #current working directory
   builddir=('dirx_%d_%s' % (case,jobname))
   log_path=('../log_%d_%s.txt' % (case,jobname))
   cwd = os.getcwd()  #current working directory
   if (stage < 2) :
      os.system('rm -rf '+builddir) #destroy directory from previous run
      os.system('mkdir '+builddir )
      os.chdir(os.path.join(cwd,builddir))
      all_opts={}
      all_opts.update(base_opts)
      all_opts.update(build_opts) # build_opts will overwrite base_opts
      all_opts.update(force_opts) # force_opts will overwrite base_opts/build_opts
      allopt_str=' '
      for k,v in all_opts.items() :
         allopt_str=allopt_str+' -D'+k+'='+v+' ' 

      logf = open(log_path, 'w')
      print( '===> OPTIONS:' )
      print( allopt_str )
      logf.write('===> OPTIONS: \n %s \n' % (allopt_str) )
      cmdconfig = 'cmake ../../.'+allopt_str
      print( '  ' )
      print( '===> configure...' )
      os.system(cmdconfig+' 2>std.err 1>std.out')
      stdout=open('std.out','r').readlines()
      stderr=open('std.err','r').readlines()

      success= True
      for line in stderr :
         if 'ERROR' in line :
            success= False
            sys.stdout.write('%s' % line)
      
      logf.write('CONFIG: STDOUT FILE:\n')
      for line in stdout :
         logf.write(line)
      
      logf.write('CONFIG: STDERR FILE:\n')
      for line in stderr :
         logf.write(line)
      
      logf.close()
      if (ntail > 0) :
         for line in stdout[-ntail:] :
            sys.stdout.write('%s' % line)
      if(not success) :
         print( '================================================ ' )
         print( ' !!!! ERROR IN CMAKE, NOT FINISHED CORRECTLY!!!! ' )
         print( '================================================ ' )
         print( '  ' )
         err.extend('caseID: %6d ,jobname: %s <=Cmake ' % (case,jobname))
         os.chdir(cwd)
         return success #=False
      else : 
         print( '================================================ ' )
         print( ' Cmake finished successfully.                    ' )
         print( '================================================ ' )
         print( '  ' )
      #endif (not success)

      #MAKE
      print( '  ' )
      print( '===> make...' )
      cmdmake = ('make -j %d VERBOSE=1' % (mpi_procs) )
      
      os.system(cmdmake+' 2>std.err 1>std.out')
      stdout=open('std.out','r').readlines()
      stderr=open('std.err','r').readlines()
      
      logf = open(log_path, 'a') #append
      
      success= False
      logf.write('\n ====> MAKE: STDOUT FILE:\n \n')
      for line in stdout :
         logf.write(line)
         if 'SUCCESS: FLUXO BUILD COMPLETE!' in line :
            success= True
      
      logf.write('\n ====> MAKE: STDERR FILE:\n \n')
      for line in stderr :
         logf.write(line)
      
      logf.close()
      if (ntail > 0) :
         for line in stdout[-ntail:] :
            sys.stdout.write('%s' % line)
      
      if (not success) :
         print( '==================================================== ' )
         print( '!!!! PROBLEM WITH BUILD, NOT FINISHED CORRECTLY!!!!! ' )
         err.extend(['caseID: %6d ,jobname: %s <=Make' % (case,jobname)])
         print( '==================================================== ' )
         print( '  ' )
         os.chdir(cwd)
         return success #=False
      else :
         print( '==================================================== ' )
         print( ' Build finished sucessfully.                         ' )                  
         print( '==================================================== ' )
         print( '  ' )
      #endif (not success)

   #endif (stage <2)
   # ...
   # RUN TEST: stage =1 and 2
   if(stage>0):

      path_to_exec=os.path.join(cwd,builddir,'bin','fluxo')
      if (not os.path.isfile(path_to_exec)) :
        err.extend(['caseID: %6d ,jobname: %s <=RUN TESTS: did not find "bin/fluxo" in builddir "%s"' % (case,jobname,builddir)])
        os.chdir(cwd)
        return False
      if(run_opts):
         for r_name,r_opts in run_opts.items():
            print('... RUN TEST "%s" ' % r_name)
            #check mandatory keys
                    
            os.chdir(os.path.join(cwd,*r_name.split('/')))

            projectname = read_prm('parameter.ini','ProjectName')
            projectnamex=('%d_%s_%s' % (case,jobname,projectname))


            try :
               # execute the code and parse the output
               success = execute(path_to_exec,'parameter.ini',projectnamex,mpi_procs=mpi_procs,ntail=ntail)
            except :
               err.extend(['caseID: %6d ,jobname: %s <=RUN %s problem in execute' % (case,jobname,r_name)])
               os.chdir(cwd)
               return False

            # RUN CHECKS
            if('test_opts' in r_opts):
               for t_name,t_opts in r_opts['test_opts'].items():
                   if 'func' in t_opts:
                      if 'f_kwargs' in t_opts:
                         [chk,msg]= t_opts['func'](**t_opts['f_kwargs'])
                      else:
                         [chk,msg]= t_opts['func']()
                      if(not chk):
                        print('   ... TEST "%s" FAILED!!!!' % t_name)
                        success=False
                        err.extend(['caseID: %6d ,jobname: %s <=RUN %s, TEST: "%s" : %s' % (case,jobname,r_name,t_name,msg)])
                      else:
                        print('   ... TEST "%s" SUCCESSFULL ...' % t_name)
            else:
               print( ' (no functions to check result specified...) ' )
            
      else:
         print( ' (no test specified...) ' )
         success = True #ignore this test
   os.chdir(cwd)
   return success

########################################################################################################
# parse a str with comma separated ranges: 1,5,10-12
def parse_range(astr):
    result = set()
    for part in astr.split(','):
       x = part.split('-')
       result.update(range(int(x[0]), int(x[-1]) + 1))
    return sorted(result)

########################################################################################################

# MAIN PROGRAM

########################################################################################################
import argparse
from job_check import job_init


assert (sys.version_info.major == 3),'python>=3.0 is needed!'

parser = argparse.ArgumentParser(description='Tool to build/compile fluxo in different configurations',\
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-p','--procs',type=int,default=1, help='    number of processors used for the make / execute (DEFAULT=1)')

parser.add_argument('-withmpi', type=int, default=1,   help='1 : DEFAULT ,compile with mpi\n'
                                                            '0 : compile without mpi  ')

parser.add_argument('-buildhdf5', type=int, default=0, help='0 : DEFAULT, use external hdf5 (modules),\n'
                                                            '1 : build hdf5 locally')

parser.add_argument('-hostname', type=str, default='', help='    cmake hostname, only needed if compiling on a cluster' )

parser.add_argument('-build_type', type=str, default='Debug', help='force a build type for all builds:\n'
                                                                    '"Debug": DEFAULT, or "Release" ' )

parser.add_argument('-case', type=str, default='0',    help='0 : DEFAULT, run all cases,\n' 
                                                            '101,102-104 : list of specific cases to run (without spaces!) ')
parser.add_argument('-tags', type=str, default='',    help='no tags: DEFAULT, run all cases,\n' 
                                                           'mhd,navierstokes: list of tags, runs only cases that have at least one matching tag')

parser.add_argument('-stage',type=int, default='0',    help='0 : DEFAULT, only build code\n'
                                                            '1 : build code and run with executable and parameterfile\n' 
                                                            '2 : only run examples (checks if code builds exist) ')

parser.add_argument('-ntail', type=int, default=5,     help='    number of last line output of cmake/make screenlog (DEFAULT=5)')




args = parser.parse_args()

runtags=[]
if(args.tags != ''):
  runtags.extend(args.tags.split(','))
cases=[]
if(args.case !='0'):
  cases = parse_range(args.case)

jobs = job_init()

# THESE OPTIONS ARE ALWAYS ENFORCED
globopts={}
globopts['FLUXO_BUILD_MPI']= ('OFF' if (args.withmpi == 0) else 'ON')
globopts['FLUXO_BUILD_HDF5']= ('OFF' if (args.buildhdf5 == 0) else 'ON')

if(len(args.hostname) > 1 ) :
   globopts['CMAKE_HOSTNAME']= args.hostname

globopts['CMAKE_BUILD_TYPE']= args.build_type

errmsg=[]

pass_args={'stage': args.stage,
           'ntail': args.ntail,
           'mpi_procs': args.procs,
          }



job_runs=[]
for j_name,job in jobs.items(): 
  tagmatch=True
  if(len(runtags)>0):
     tagmatch = any([(tag in job['tags']) for tag in runtags])
  casematch=True
  if(len(cases)>0):
     casematch= (job['case'] in cases)
  if(casematch and tagmatch):
     job_runs.extend([j_name])
     print('\n \n caseID: %s jobname: %s, tags: %s \n '% (job['case'],j_name,job['tags']))
     stat = test_fluxo(jobname=j_name,force_opts=globopts,err=errmsg,**pass_args,**job)

#### GLOBAL CHECK

if(len(job_runs) == 0):
   print( ' NO MATCHING JOBS FOUND FULFILLING -case AND -tags ARGUMENTS!' )
   sys.exit(100)

if(len(errmsg) > 0 ) :
   print( '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' )
   print( '!!!!!!!    WARNING, following jobs failed:     !!!!!!!!' )
   print( ' ' )
   for line in errmsg :
      print( '--> '+line )
   print( ' ' )
   print( '... see log.caseID_project files' )
   print( '   and dirx_caseID_project folders.' )
   print( ' rerun this script with -case option to run again spefic caseIDs.' )
   print( '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' )
   sys.exit(100)
else :
   print( '/////////////////////////////////////////////////////////' )
   print( ' ' )
   if ( args.stage == 0) :
     print( ' ==> ALL BUILDS SUCCESSFULL!' )
   if ( args.stage == 1) :
     print( ' ==> ALL BUILDS AND RUNS SUCCESSFULL!' )
   if ( args.stage == 2) :
     print( ' ==> ALL RUNS SUCCESSFULL!' )
   print( ' ' )
   print( '/////////////////////////////////////////////////////////' )
   sys.exit(0)
