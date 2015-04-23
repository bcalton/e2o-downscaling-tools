@echo off
rem normally run a whole year, for this example only 10 days!
rem python ..\..\e2o_dstools\e2o_radiation.py -S 1 -E 366 D highResDEM\wflow_dem.map -d lowResDEM\dem.tif -O ouput_rad -s 5 -e 22

rem only use fast (but lest accurate 180 min internal timesteps)
python ..\..\e2o_dstools\e2o_radiation.py -S 1 -E 5 -D highResDEM\wflow_dem.map -d lowResDEM\dem.tif -O ouput_rad -s 5 -e 22 -T180 -l DEBUG

