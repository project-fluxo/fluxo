import os
import sys
from job_def import job_definition
#import helpers


###########################################################################
# run job_definition and check all the jobs!
###########################################################################
def job_init():
   
   jobs=job_definition()

   jobs_correct=True
   for j_name,job in jobs.items():
      errmsg=[]
      #check mandatory keys in job
      for mand_key in ['case','tags','build_opts']:
         if (not (mand_key in job)) : 
            errmsg.append('mandatory "%s" key missing in job' % (mand_key))
      if('run_opts' in job):
         for r_name,r_opts in job['run_opts'].items():
            for mand_key in ['tags']:
               if (not (mand_key in r_opts)) : 
                  errmsg.append('mandatory "%s" key missing in run_opts of %s' % (mand_key,r_name))
            #check mandatory keys in run_opts
            rdir=os.path.join(*r_name.split('/'))
            if(not os.path.isdir(rdir)):
               errmsg.append('"%s" directory does not exist ' % (rdir))
            else:
               rprm_file=os.path.join(rdir,'parameter.ini')
               if(not os.path.isfile(rprm_file)):
                  errmsg.append('parameter.ini does not exist in directory "%s" ' % (rdir))
   
            if('test_opts' in r_opts):
               for t_name,t_opts in r_opts['test_opts'].items():
                  #check mandatory keys in test_opts
                  if (not ('func' in t_opts)) : 
                     errmsg.append('"%s" key missing in job["run_opts"]["test_opts"], name "%s" ' % (mand_key,t_name))
   
      if(len(errmsg) >0):
         jobs_correct=False
         print('!!!! PROBLEM WITH  DEFINITION OF JOB "%s":'%(j_name))
         for line in errmsg :
            print( "!!!!===> "+line )
         print(job)
         
   if(jobs_correct):
      print("===>  ALL JOBS CHECKED")
      return jobs
   else:
      sys.exit(100)


##test
#jobs=job_init()
#ff= jobs['navierstokes_type1_br1_GL']['run_opts']['runs/navst/freestream/conforming']['test_opts']['err_Linf']
#
#[stat,msg]=ff['func'](**ff['f_kwargs'])
#print(stat,msg)
