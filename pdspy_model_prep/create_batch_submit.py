import os

def create_batch_submit(source_name, line_name, source_dir, disk_type, ncpu):

    model_path = source_dir + disk_type + '_' + line_name + '/'
    
    with open(model_path + 'runmodel.csh', 'w') as runmodel:
        run_command = "xvfb-run mpiexec -np {0} flared_model_nested.py -o \
{1} -n 1 --ftcode galario-unstructured".format(ncpu, source_name)
        runmodel.write(run_command)
