#!/usr/bin/env ipython3
import socket
import getpass

def get_basepath():
    host_name = socket.gethostname()
    user_name = getpass.getuser()
    basepath = '/home/psteger/sci/darcoda/gravimage/'
    if 'darkside' in host_name:
        basepath = '/home/ast/read/dark/darcoda/gravimage/'
    elif ('lisa' in host_name) and ('hsilverw' in user_name):
        basepath = '/home/hsilverw/LoDaM/darcoda/gravimage/'
    elif ('lisa' in host_name) and ('sofia' in user_name):
        basepath = '/home/sofia/darcoda/gravimage/'
    return basepath
## \fn get_basepath()
# return basepath, depending on machine

def get_machine():
    host_name = socket.gethostname()
    user_name = getpass.getuser()
    if 'darkside' in host_name:
        machine = 'darkside'
    elif 'science' in host_name:
        machine = 'science'
    elif ('lisa' in host_name) and ('login' in host_name) and ('hsilverw' in user_name):
        machine = 'lisa_HS_login'
    elif ('lisa' in host_name) and ('login' not in host_name) and ('hsilverw' in user_name):
        machine = 'lisa_HS_batch'
    elif ('lisa' in host_name) and ('login' in host_name) and ('sofia' in user_name):
        machine = 'lisa_SS_login'
    elif ('lisa' in host_name) and ('login' not in host_name) and ('sofia' in user_name):
        machine = 'lisa_SS_batch'
    return machine
## \fn get_machine()
# return string for machine


# HS:
# check whether following code pasted from import_path.py.set_geometry() needs be factored in

# if machine == 'darkside':
#     basepath = '/home/ast/read/dark/darcoda/gravimage/programs/'
# elif machine == 'lisa_HS_login':
#     basepath = '/home/hsilverw/LoDaM/darcoda/gravimage/programs/'
# elif machine == 'lisa_SS_login':
#     basepath = '/home/sofia/darcoda/gravimage/programs/'
# elif machine == 'lisa_HS_batch' or machine == 'lisa_SS_batch':
#     scratch_space = os.getenv("TMPDIR")
#     basepath = scratch_space + '/darcoda/gravimage/programs/'
#print('basepath = ', basepath)
