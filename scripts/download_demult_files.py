'''
Author: Mudra Hegde
Email: mhegde@broadinstitute.org
'''

import argparse, os, sys


def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('--username', type=str, help='Username provided by GP')
	parser.add_argument('--password', type=str, help='Password provided by GP')
	parser.add_argument('--url', type=str, help='URL provided by GP')
	parser.add_argument('--output-folder-name', type=str, help='Name of output folder to save demultiplexed files in')
	parser.add_argument('--output-folder-path', type=str, help='Path to output folder to save demultiplexed files in'
															   '(Do not include the name of the folder itself)')
	parser.add_argument('--make-backup', type=str, default='N', help='Y/N to make backup of files on server; Default:N')
	return parser


def get_inputs():
	args = get_parser().parse_args()
	username = args.username
	password = args.password
	url = args.url
	outputname = args.output_folder_name
	outputpath = args.output_folder_path
	make_backup = args.make_backup
	return username, password, url, outputname, outputpath, make_backup


def download_locally():
	username, password, url, outputname, outputpath, make_backup = get_inputs()
	dir_cmd = 'mkdir '+outputpath+'/'+outputname
	os.system(dir_cmd)
	cmd = 'wget --tries=10 --continue --mirror --user '+ username +' --password '+ password +\
		  ' --no-check-certificate '+ url +' -P '+outputpath+'/'+outputname
	try:
		sys.path.append('/usr/local/bin/')
		os.system(cmd)
	except:
		sys.exit()
	if make_backup.upper() == 'Y':
		print("Making a backup of the files on the server...")
		make_stable_backup(outputname, outputpath)
	return


def make_stable_backup(outputname, outputpath):
	path = '/Volumes/gpp_share/old/Mhegde/DemultiplexedData/Demultiplexed/'
	dir_cmd = 'mkdir '+path+'/'+outputname
	os.system(dir_cmd)
	copy_cmd = 'cp -r '+outputpath+'/'+outputname+' '+path+'/'+outputname
	os.system(copy_cmd)
	return


if __name__ == '__main__':
	download_locally()