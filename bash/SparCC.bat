echo OFF

:: Define pathways for input files, outputfiles and to python
set SparCCPackage=C:\Users\Jule\anaconda3\pkgs\sparcc-0.1.0-0\python-scripts\
set InputPath=C:\Users\Jule\sciebo\AG_Bonkowski\Projects\ContrastingActivityPatterns\DataToAnalyse\NetworkData\Raw\
set OutputPath=C:\Users\Jule\sciebo\AG_Bonkowski\Projects\ContrastingActivityPatterns\DataToAnalyse\NetworkData\SparCC\

:: Create directory for output files
if not exist "%OutputPath%" (
   mkdir "%OutputPath%"
)

:: Define Data type and WWTP compartments
set data_type=DNA,RNA
set compartments=INF,DNF,NFC,EFF

:: Loop iterates over the data types (DNA/RNA)
for %%d in (%data_type:,= %) do (
	echo data type: %%d
	
	:: Loop iterates over the WWTP compartments
	for %%c in (%compartments:,= %) do (
	
		echo compartment: %%c

		echo Step 1: Compute correlations
		
		python %SparCCPackage%SparCC.py %InputPath%%%d_%%c.txt -i 20 -c %OutputPath%%%d_%%c_SparCC.txt > %OutputPath%%%d_%%c_Logfile_corr.log

		echo Step 2: Re-sample data
		
		if not exist %OutputPath%Resampling\%%d_%%c (
			mkdir %OutputPath%Resampling\%%d_%%c
		)
		
		python %SparCCPackage%MakeBootstraps.py %InputPath%%%d_%%c.txt -n 100 -t permutation_#.txt -p %OutputPath%Resampling\%%d_%%c\
		
		echo Step 3: Compute correlations with re-sampled data

		if not exist %OutputPath%Bootstraps\%%d_%%c (
			mkdir %OutputPath%Bootstraps\%%d_%%c
		)

		for /l %%a in (0,1,99) do (
			python %SparCCPackage%SparCC.py %OutputPath%Resampling\%%d_%%c\permutation_%%a.txt -c %OutputPath%Bootstraps\%%d_%%c\bootstraps_%%a.txt >> %OutputPath%%%d_%%c_Logfile_bootstraps.log
		)

		echo Step 4: Compute pseudo p-values
		python %SparCCPackage%PseudoPvals.py %OutputPath%%%d_%%c_SparCC.txt %OutputPath%Bootstraps\%%d_%%c\bootstraps_#.txt 100 -o %OutputPath%%%d_%%c_pvals_two_sided.txt -t two_sided >> %OutputPath%%%d_%%c_Logfile_pvals.log

	)
)
