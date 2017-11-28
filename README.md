# README #

* sqa is developed by James Kneller. This version is edited by Yonglin Zhu for neutrino oscillation calculation with matter and neutrino density profile from simualtions. Sherwood and Yonglin Zhu implemented some other physics (scattering, etc.) on the fly. 

### What is this repository for? ###


### How do I get set up? ###

* sqa need to be compiled with mstl library()
1. Get mstl ready:
- On HPC@NC State: Check if you have access to /home/jpknelle/mstl/, if you do, then you can keep using the makefiles.
- On /physics: check if you have access to /physics/jpknelle/mstl/, if you do, then you can keep using the makefile
s.
- Some where else: mstl need to be compiled again if it is used somewhere else. (To compile mstl, in mstl/lib, edit path in *txt and run libraries.sh)
2. Get sqa ready: 
- With icc, if you work on HPC or /physics, just compile with the makefiles
- With gcc, if you work on HPC or /physics, just compile with the makefiles
3. Run sqa
- On HPC@NC State: find the sumbit file, which should look like:
~~~~~~~~~~~~~~~~~~
#!/bin/csh
#BSUB -R rusage[mem=2048] 
#BSUB -o ./outerr/%J.out 
#BSUB -e ./outerr/%J.err
#BSUB -W 50 
#BSUB -q nucastro
#BSUB -n 4 
#BSUB -R span[hosts=1]

./sqa2.psma.x sqa2.psma.run

~~~~~~~~~~~~~~~~~
- -R is the ram usage; -o/-e is the out and err files generate by HPC;
- -W set the upper limit of your running time, it may change the time you need to wait in line.
- -q is the queue name, nucastro is our queue with 128 cores. Most of them should be free.
- -n is the number of cores you request. I am not sure the upper limit. You may ask Jim about that
- -R how you want the cores distribute, whether on one host or multi-host
- Before you start the first run, you need to change output directory in .run file. You should be able to 
- make a directory of your own in /gpfs_common/share01/jpknelle/, that is where I put the output. However, it will be cleaned every two weeks, if you do not touch that directory. 
- Once you get your sumbit file and .run ready, type:
~~~~~~~~
bsub < submit
~~~~~~~~
- You will be assign a job id. To check it, `bjobs`.

- On reines/cowan: Similiar to HPC, but you do not need the sumbit file, just run it with `./sqa2.psma.x sqa2.psma.run`

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact
