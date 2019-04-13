In ~/.ssh/config, put:
Host thor
User 4100484
HostName thor.science.uu.nl

In the cluster, in your user folder, make sure there exists a file .sge_request, containing only "-cwd", without quotation marks.
This makes it so that when submitting a job, the output data goes to the Current Working Directory.
Otherwise manually adding -cwd to each job becomes a chore, and not doing this clutters up your user folder.

Note: don't use fflush on programs in the cluster. There is only limited bandwidth for communication between nodes and the headnode.
Also don't output too much data regardless, and make sure to remove data you don't need anymore.

Commands:

Connect to the cluster:
ssh -X thor
(-X to be able to use programs which pop out of the terminal)

Show programs I run:
qstat

Show what everyone is running:
qstat -u "*"

Show data per machine:
qstat -f

Submit a job to the cluster:
qsub filename.bat
(use .bat extension to prevent problems)
filename.bat then contains:
cd directory
./executable.exe > out.dat

Stop a job:
qdel job_id
(use qstat to find the job_id)

Copy data from the cluster:
Go to a terminal on your own machine, then you have two options.
1)
scp thor:~/out.dat .
(copies ~out.dat on thor to your local pwd (.))

2)
sftp thor
(now you can browse on the cluster. (cd, ls, pwd and so on work, albeit slowly))
get out.dat
(copies out.dat from the cluster's pwd to your terminal's pwd)
put script.sh
(copies script.sh from your terminal's pwd to the clusters pwd)
lls, lcd, lpwd
(local ls, cd, pwd. i.e. do these commands on your local machine)
