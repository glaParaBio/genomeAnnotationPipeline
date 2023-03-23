rule install_augustus:
    output:
        os.path.join(INSTALL_PATH, 'augustus/bin/augustus'),
        directory(os.path.join(INSTALL_PATH, 'augustus/config')),
    params:
        git_commit= 'd1fd097',
        INSTALL_PATH= INSTALL_PATH,
    shell:
        r"""
        rm -rf Augustus
        git clone https://github.com/Gaius-Augustus/Augustus.git
        cd Augustus
        git checkout {params.git_commit}

        git checkout common.mk # Reset to original before editing it
        > common.tmp
        echo "
        INCLUDE_PATH_ZLIB        := -I{params.INSTALL_PATH}/include
        INCLUDE_PATH_BOOST       := -I{params.INSTALL_PATH}/include
        INCLUDE_PATH_LPSOLVE     := -I{params.INSTALL_PATH}/include/lpsolve
        INCLUDE_PATH_SUITESPARSE := -I{params.INSTALL_PATH}/include/suitesparse/
        INCLUDE_PATH_GSL         := -I{params.INSTALL_PATH}/include
        INCLUDE_PATH_SQLITE      := -I{params.INSTALL_PATH}/include
        INCLUDE_PATH_BAMTOOLS    := -I{params.INSTALL_PATH}/include/bamtools
        INCLUDE_PATH_HTSLIB      := -I{params.INSTALL_PATH}/include/htslib
        INCLUDE_PATH_SEQLIB      := -I {params.INSTALL_PATH}/include/SeqLib -I{params.INSTALL_PATH}/include/htslib -I{params.INSTALL_PATH}/include/jsoncpp
        LIBRARY_PATH_LPSOLVE     := -L{params.INSTALL_PATH}/lib/lp_solve/ -Wl,-rpath,{params.INSTALL_PATH}/lib/lp_solve/

        LIBRARY_PATH_ZLIB        := -L{params.INSTALL_PATH}/lib -Wl,-rpath,{params.INSTALL_PATH}/lib
        LIBRARY_PATH_BOOST       := -L{params.INSTALL_PATH}/lib -Wl,-rpath,{params.INSTALL_PATH}/lib
        LIBRARY_PATH_SUITESPARSE := -L{params.INSTALL_PATH}/lib -Wl,-rpath,{params.INSTALL_PATH}/lib
        LIBRARY_PATH_GSL         := -L{params.INSTALL_PATH}/lib -Wl,-rpath,{params.INSTALL_PATH}/lib
        LIBRARY_PATH_SQLITE      := -L{params.INSTALL_PATH}/lib -Wl,-rpath,{params.INSTALL_PATH}/lib
        LIBRARY_PATH_BAMTOOLS    := -L{params.INSTALL_PATH}/lib -Wl,-rpath,{params.INSTALL_PATH}/lib
        LIBRARY_PATH_HTSLIB      := -L{params.INSTALL_PATH}/lib -Wl,-rpath,{params.INSTALL_PATH}/lib
        LIBRARY_PATH_SEQLIB      := -L{params.INSTALL_PATH}/lib -Wl,-rpath,{params.INSTALL_PATH}/lib
        " >> common.tmp
        cat common.mk >> common.tmp
        mv common.tmp common.mk

        make -j 8 MYSQL=false augustus auxprogs

        mkdir -p {params.INSTALL_PATH}/augustus/bin
        mv bin/* {params.INSTALL_PATH}/augustus/bin/
        cp -r scripts/* {params.INSTALL_PATH}/augustus/bin/
        cp -r config {params.INSTALL_PATH}/augustus/

        # Some basic tests:
        {params.INSTALL_PATH}/augustus/bin/augustus --version
        {params.INSTALL_PATH}/augustus/bin/autoAug.pl --help
        """


rule install_braker:
    # Braker is on bioconda but it may be outdated so we install it ourselves
    output:
        os.path.join(INSTALL_PATH, 'bin/braker.pl'),
    params:
        git_commit= '1af9bdc',
        INSTALL_PATH= INSTALL_PATH,
    shell:
        r"""
        rm -rf BRAKER
        git clone https://github.com/Gaius-Augustus/BRAKER
        cd BRAKER
        git checkout {params.git_commit}

        chmod +x scripts/*.pl scripts/*.py
        cp scripts/*.pl {params.INSTALL_PATH}/bin/
        cp scripts/*.pm {params.INSTALL_PATH}/bin/
        cp scripts/*.py {params.INSTALL_PATH}/bin/
        cp -r scripts/cfg/ {params.INSTALL_PATH}/bin/

        rm -rf BRAKER
        braker.pl --help > /dev/null
        """


rule install_genemark:
    input:
        gm= config['genemark_tar_gz'],
    output:
        touch('genemark_installation.done'),
    params:
        gmpath= GENEMARK_PATH,
    shell:
        r"""
        gm=`basename {input.gm} .tar.gz`
        rm -rf {params.gmpath}
        cp {input.gm} ${{gm}}.tar.gz
        tar zxf ${{gm}}.tar.gz
        mv ${{gm}} {params.gmpath} 
        rm -rf ${{gm}}.tar.gz

        cd {params.gmpath}

        # Use perl on PATH
        perl change_path_in_perl_scripts.pl '/usr/bin/env perl'
        
        # Test installation
        ./check_install.bash
        """


