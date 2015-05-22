@echo off
rem normally run a whole year, for this example only 10 days!
rem python ..\..\e2o_dstools\e2o_radiation.py -S 1 -E 366 D highResDEM\wflow_dem.map -d lowResDEM\dem.tif -O output_rad -s 5 -e 22

rem only use fast (but lest accurate 180 min internal timesteps)
rem start python ..\..\e2o_dstools\e2o_radiation.py -S 1 -E 366 -D highResDEM\wflow_dem.map -d lowResDEM\dem.tif -O output_rad -s 5 -e 22 -T180 -l DEBUG
rem only use fast (but lest accurate 180 min internal timesteps)
start python ..\..\e2o_dstools\e2o_radiation.py -S 1 -E 159 -D highResDEM\wflow_dem.map -d lowResDEM\dem.tif -O output_rad_linke -s 5 -e 22 -T180 -l DEBUG -L ../../data/linke
start python ..\..\e2o_dstools\e2o_radiation.py -S 160 -E 366 -D highResDEM\wflow_dem.map -d lowResDEM\dem.tif -O output_rad_linke -s 5 -e 22 -T180 -l DEBUG -L ../../data/linke

