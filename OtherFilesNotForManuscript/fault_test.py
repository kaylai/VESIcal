import faulthandler; faulthandler.enable() #to handle segmentation faults thrown by MELTS/thermoengine
import VESIcal as v 

myfile = v.ExcelFile('Statia_MI_Cooper_et_al_5.xlsx')

satPs = myfile.calculate_saturation_pressure(temperature="Temp")
print(satPs)