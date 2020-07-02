import VESIcal as v
import anvil.server
import anvil.media

anvil.server.connect("KYKUMVNXIZ23YPRPDDGL5DM6-E4XIIJ7V23RSDNS5")

@anvil.server.callable
def import_ExcelFile(file):
	with anvil.media.TempFile(file) as filename:
		myfile = v.ExcelFile(filename)

	return(str(filename) +" has been loaded!")