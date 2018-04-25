arrayjob_for_cmds.py jobs.json "P12_wcsmin_autocover" 4 32GB 22:00:00 11 > fine_job.pbs
arrayjob_for_cmds.py jobs_coarser_p8.json "P8_wcsmin_autocover" 4 16GB 1:00:00 > coarser_p8_job.pbs
arrayjob_for_cmds.py jobs_coarser_p8.json "P8_wcsmin_autocover" 4 32GB 1:00:00 > coarser_p8_job_32GB.pbs
arrayjob_for_cmds.py jobs_coarser_p8_ob0.2.json "P8_wcsmin_ac_ob0.2" 4 16GB 1:00:00 > coarser_p8_ob2_job.pbs

arrayjob_for_cmds.py jobs_coarser_p8_1_12_2017.json "P8_wcsmin_autocover" 4 16GB 1:00:00 > coarser_p8_job_1_12_2017.pbs
