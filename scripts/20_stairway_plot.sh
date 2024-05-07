#!/bin/bash

cd /piec1/software/stairway-plot-v2/stairway_plot_v2.1.1/

for file in *.blueprint;
	do
		/usr/bin/java -cp /piec1/software/stairway-plot-v2/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder /piec1/software/stairway-plot-v2/stairway_plot_v2.1.1/$file

		chmod 777 /piec1/software/stairway-plot-v2/stairway_plot_v2.1.1/$file.sh

		/piec1/software/stairway-plot-v2/stairway_plot_v2.1.1/$file.sh

		chmod 777 /piec1/software/stairway-plot-v2/stairway_plot_v2.1.1/$file.plot.sh

		/piec1/software/stairway-plot-v2/stairway_plot_v2.1.1/$file.plot.sh
	done




