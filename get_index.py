def get_index(server):

    import os
    import config

    # check if local Argo index file is up to date
    if not os.path.exists(config.fidx): # download Argo index file
        cmd = 'ftp -V '+server+config.fidx+'.gz'
        print cmd
        print 'Downloading Argo index file. This may take some time... '
        os.system(cmd)  # get updated argo index file
        cmd = 'gunzip -f '+config.fidx+'.gz'
        os.system(cmd)
    else: # check the one you have is up-to-date
        cmd = 'curl -r 0-500 '+server+'ar_index_global_prof.txt > tmp_header.txt'
        print 'Checking timestamp of Argo index file on server... '
        print cmd
        os.system(cmd)

        f = open('tmp_header.txt','r')
        for line in f:
            if 'Date of update' in line: dlineServer=line
        f.close()

        f = open(config.fidx,'r')
        for line in f:
            if 'Date of update' in line: dlineLocal=line
        f.close()

        uptodate = dlineServer == dlineLocal

        cmd = 'rm tmp_header.txt'
        os.system(cmd)

        if not uptodate:
            cmd = 'ftp -V '+server+config.fidx+'.gz'
            print 'Updating local Argo index file. This may take some time... '
            print cmd
            os.system(cmd)  # get updated argo index file
            cmd = 'gunzip -f '+config.fidx+'.gz'
            os.system(cmd)
        else:
            print 'The local Argo index file is up-to-date.'
