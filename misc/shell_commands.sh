#!/usr/bin/env bash
# -*- coding: utf-8 -*-


#mahti to gtb
rsync -av -e ssh mahti:/scratch/project_2014077/geolifeclef/hmsc/ /Users/gtikhono/DATA/geolifeclef/geolifeclef-2025/hmsc/ ; rsync -av /Users/gtikhono/DATA/geolifeclef/geolifeclef-2025/hmsc/ -e ssh gtb:/home/gt/DATA/geolifeclef-2025/hmsc/

#gtb to mahti
rsync -av -e ssh gtb:/home/gt/DATA/geolifeclef-2025/hmsc/init/ /Users/gtikhono/DATA/geolifeclef/geolifeclef-2025/hmsc/init/ ; rsync -av /Users/gtikhono/DATA/geolifeclef/geolifeclef-2025/hmsc/init/ -e ssh mahti:/scratch/project_2014077/geolifeclef/hmsc/init/
rsync -av --exclude={"SatelitePatches","SateliteTimeSeries-Landsat","BioclimTimeSeries","hmsc/pred"} -e ssh gtb:/home/gt/DATA/geolifeclef-2025/  /Users/gtikhono/DATA/geolifeclef-2025/ ; rsync -av --exclude={"SatelitePatches","SateliteTimeSeries-Landsat","BioclimTimeSeries","hmsc/pred"} /Users/gtikhono/DATA/geolifeclef-2025/ -e ssh mahti:/scratch/project_2014077/geolifeclef/


