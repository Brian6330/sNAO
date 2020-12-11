## Summer North Atlantic Oscillation:
### Effects of Summer NAO on Central Europe and the North American East Coast
#### A project by Louis Zwyssig (louis.zwyssig@students.unibe.ch) & Brian Schweigler (brian.schweigler@students.unibe.ch)

### What is the North Atlantic Oscillation (NAO)?
The NAO is the only teleconnection pattern evident throughout the year in the Northern Hemisphere (NH). It is described as a redistribution of atmospheric mass between the Arctic and the subtropical Atlantic. 

The swings from one phase to another produce large changes that can be felt in the weather for the Atlantic and neighboring continents. These changes include: 
mean wind speed and direction, heat, moisture, as well as intensity and frequency of storms.

### What is the Summer North Atlantic Oscillation (sNAO)? Why is it important?
In the NH, the sNAO is defined during the high summer months of July and August, as the NAO in those months. Centers of Action are similar to NAO but should be adjusted for increased precision. 

The sNAO exerts a strong influence on northern European rainfall, temperature, and cloudiness through changes in the position of the North Atlantic storm track. It is therefore key in generating summer climate extremes (flooding, drought, and heat stress) in northwestern Europe. 

### Methodology:
Due to the  pressure difference at geographic locations, the sea-level pressure (slp) must be standardized to be compared. Thus, to define the sNAO, standardized slp was used, based on Folland et al. (2009) and their defined region of 25.8°-70.8°N and 70.8°W-50.8°E, with the leading principal component (pc) of slp anomalies for high summer months July and August.

The centers of action found through the leading pc were used as basis for regression and correlation of precipitation, temperature, and geopotential height at 500 hPa over the region of 22°-82°N, 78°W-52°E.

These values were then plotted to facilitate the analysis of possible effects of the sNAO on our weather.
All scripts and plots were created with R (found on: https://github.com/Brian6330/sNAO)

### Limitations:
Data in EKF400 can have random and systematic errors; “Observations [used in EKF400] also have errors, e.g., due to instrument changes, reporting errors, or non-climatic signals in proxies.”
- By using a range of 400 years, the effect of random errors is minimized, but not for systematic errors, which may be even more prevalent in early data. 

The impact of other large-scale atmospheric processes is unclear. 
- We focus on the sNAO to ascertain its impact, but it is impossible to observe a process disparate from others (e.g., sNAO is influenced by ENSO (TODO Louis: insert reference))

Falsely restraining spatial effects of teleconnected processes by limiting on small geographic region.
- As the aim is to find impacts on Central European and East-North American weather, it is valid to concentrate solely on those regions, despite the global resolution of the dataset.

### Noticeable Effects of the sNAO index on central European weather:
Positive  sNAO Phase: 	below average Precipitation, 	above average Temperature
Negative sNAO Phase: 	above average Precipitation, 	below average Temperature 

### Findings 
- There is a clear relationship between the strength of sNAO phenomenon and European and North-East American precipitation and temperature, which is consistent with findings of Folland et al.
- Center of action over the Mediterranean stronger in EKF400 dataset than Folland et al. 
- The results are supported by case study of London temperature and precipitation anomalies in extreme sNAO summers (e.g. 1783, 1955).

_We would like to thank Prof. Dr. Stefan Brönnimann for his assistance with the R-scripts as well access to the EKF400 dataset._

### References:
1. Barnston, A. G., Livezey, R. E. 1987. Classification, Seasonality and Persistence of Low-Frequency Atmospheric Circulation Patterns. Monthly Weather Review. 115(6). 1083-1126.
2. Folland, C. K., Knight, J., Linderholm, H. W., Fereday, D., Ineson, S., Hurrell, J. W. 2009. The Summer North Atlantic Oscillation: Past, Present, and Future. Journal of Climate. 22(5). 1082-1103.
3. Jones, P. D., Jonsson, T., Wheeler, D. 1997. Extension to the North Atlantic Oscillation Using Early Instrumental Pressure Observations from Gibraltar and South-West Iceland. International Journal of Climatology. 17. 1433-1450.
4. Franke, J., Brönnimann, S., Bhend, J., Brugnara, Y. 2017. A monthly global paleo-reanalysis of the atmosphere from 1600 to 2005 for studying past climatic variations. Scientific Data. (4)170076. 1-19.
