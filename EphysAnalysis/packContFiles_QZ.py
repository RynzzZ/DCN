import sys
sys.path.insert(0, 'D:\\Cerebellum_DCN_Project\\Github\\analysis-tools\\Python3')
from OpenEphys import pack_2



# get input arguments
file_dir = 'Z:\\Qianyun\\DCN\\Data\\20210929_000\\2021-09-29_18-00-22\\Record Node 101'
source = 100
fs = 30000
highpass = 300
dref = 'ave'

connected_channels = 'all'

# format input arguments appropriately
fs = int(fs)
highpass = int(highpass)
if dref=='None': dref=None


print('running pack_2...')
# WITH highpass
pack_2(folderpath = file_dir, filename = '', channels = 'all', chprefix = '', highpass=highpass, fs = fs,
       dref = dref, session = '0', source = '100')
# WITHOUT highpass
# pack_2(folderpath = file_dir, filename = '', channels = 'all', chprefix = '', fs = fs, dref = dref, session = '0', source = '100')