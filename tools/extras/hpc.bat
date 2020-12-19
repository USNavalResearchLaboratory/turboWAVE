:: provide a common interface for remote HPC machines (DOS version)
:: WARNING: assumes PuTTY is installed (could use openssh, but HPC centers like PuTTY)
:: Examples:
:: 1) hpc cori login // login to Cori
:: 2) hpc cori th // sftp session starting in Cori home directory
:: 3) hpc cori t // sftp session starting in Cori scratch directory
:: 4) hpc onyx get "*.inp" // get all .inp files from Onyx scratch directory
:: 5) hpc onyx geth "*.dat" // get all .dat files from Onyx home directory
:: 6) hpc onyx put "*.inp" // put all .inp files in Onyx scratch directory
:: 7) hpc onyx puth "*.inp" // put all .inp files in Onyx home directory
:: 8) hpc cori update // put all files in working directory into ~/turbowave
:: note use of quotes enables wildcards

@ECHO OFF

IF "%1"=="cori" (
		SET lopt=
		SET login=gordond@cori.nersc.gov
		SET workdir=/global/cscratch1/sd/gordond
		SET homedir=/global/homes/g/gordond
)
IF "%1"=="onyx" (
		SET lopt=
		SET login=gordon@onyx.erdc.hpc.mil
		SET workdir=/p/work/gordon
		SET homedir=/p/home/gordon
)
IF "%1"=="must" (
		SET lopt=
		SET login=gordon@mustang.afrl.hpc.mil
		SET workdir=/p/work1/gordon
		SET homedir=/p/home/gordon
)
IF "%1"=="thun" (
		SET lopt=
		SET login=gordon@thunder.afrl.hpc.mil
		SET workdir=/workspace/gordon
		SET homedir=/home/gordon
)
IF "%1"=="water" (
		SET lopt=
		SET login=gordon@ppdwater.nrl.navy.mil
		SET workdir=/home/gordon/Run
		SET homedir=/home/gordon
)

IF "%login%"=="" (
	ECHO Unknown machine name.
)

IF NOT "%2"=="login" (
	IF NOT "%2"=="put" (
		IF NOT "%2"=="get" (
			IF NOT "%2"=="puth" (
				IF NOT "%2"=="geth" (
					IF NOT "%2"=="update" (
						IF NOT "%2"=="t" (
							IF NOT "%2"=="th" (
								ECHO Invalid action
								ECHO Supported actions: login,put,puth,get,geth,update,t,th
							)
						)
					)
				)
			)
		)
	)
)

IF "%2"=="login" (
	ECHO Sorry we don't know how, please run PuTTY from the start menu.
)

IF "%2"=="home" (
	psftp %login%
)

IF "%2"=="work" (
	psftp %login%
)

IF "%2"=="put" (
		IF "%3"=="" (
			ECHO Missing argument to put
		) ELSE (
			pscp %3 %login%:%workdir%
		)
)

IF "%2"=="puth" (
		IF "%3"=="" (
			ECHO Missing argument to puth
		) ELSE (
			pscp %3 %login%:%homedir%
		)
)

IF "%2"=="get" (
		IF "%3"=="" (
			ECHO Missing argument to get
		) ELSE (
			pscp %login%:%workdir%/%3 .
		)
)

IF "%2"=="geth" (
		IF "%3"=="" (
			ECHO Missing argument to geth
		) ELSE (
			pscp %login%:%homedir%/%3 .
		)
)

IF "%2"=="t" (
		ECHO Sorry we don't know how to set remote directory to other than home.
		ECHO You will have to use th.
)

IF "%2"=="th" (
		psftp %login%
)

IF "%2"=="update" (
		pscp * %login%:%homedir%/turbowave/
)
