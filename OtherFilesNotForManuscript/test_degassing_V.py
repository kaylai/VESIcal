import VESIcal as v 
import random

myfile = v.ExcelFile('testDataSets/example_data.xlsx')
samp = myfile.get_sample_oxide_comp('10*')

bulk_comp = samp 
temperature = 1200
temp=temperature
pressure = 'saturation'

print("Fractionate vapor = 0")
print(v.calculate_degassing_path(sample=samp, temperature=temperature, fractionate_vapor=0, init_vapor=0).result)
print("\n")

while True:
	fractionate_vapor=random.uniform(0, 1)
	print("Fractionate vapor = " + str(fractionate_vapor))
	init_vapor = 0.0

	"""Calculate open, closed, and closed + 2 wt% initial vapor"""
	#closed_df = v.calculate_degassing_path(sample=samp, temperature=temp).result
	#open_df = v.calculate_degassing_path(sample=bulk_comp, temperature=temp, fractionate_vapor=1.0).result
	#half_df = v.calculate_degassing_path(sample=bulk_comp, temperature=temp, fractionate_vapor=fractionate_vapor, init_vapor=init_vapor).result
	#exsolved_df = v.calculate_degassing_path(sample=bulk_comp, temperature=temp, init_vapor=2.0).result

	"""Calculate closed-system degassing starting from a pressure of 2000 bars"""
	#start2000_df = v.calculate_degassing_path(sample=bulk_comp, temperature=temp, pressure=2000.0).result

	print(v.calculate_degassing_path(sample=samp, temperature=temperature, fractionate_vapor=fractionate_vapor, init_vapor=init_vapor).result)