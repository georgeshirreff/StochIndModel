def readParam(param_name,param_file,set,out):
	f=open(param_file,"r")
	param_value=float('NaN')
	while 1:
		try:
			line = f.next()
			if line.find(param_name) != -1:
				param_values=line.split("\t")
				param_value=float(param_values[set])
				break
		except StopIteration:
			break
	
	out.write(param_name+"\t"+str(param_value)+"\n")
	
	return param_value
	
def readString(param_name,param_file,set):
	f=open(param_file,"r")
	param_value="NA"
	while 1:
		try:
			line = f.next()
			if line.find(param_name) != -1:
				param_values=line.split("\t")
				param_value=param_values[set]
				param_value=param_value.rstrip()
				break
		except StopIteration:
			break
	
	print(param_name+"\t"+param_value)
	
	return param_value	
	
