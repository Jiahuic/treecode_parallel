This program is built on mvapich2/2.0-gcc-4.9.1
To run the program:
when running jobs on a dedicated parallel computer (like a laptop, workstation or dedicated cluster), parallel jobs and processes are not regulated through a queueing system.
The calling syntax of mpiexec is:

{% highlight bash %}
mpiexec [mpiexec_options] program_name [program_options]
{% endhighlight %}

submit.job is a batch submission script for SLURM queueing system



