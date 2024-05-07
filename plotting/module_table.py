import sys

class TexTabler:
	def __init__(self, file_destination:str, n_colums:int):
		self.destination = file_destination
		self.n_colums = n_colums
		self.text = ""


	def add_row(self, cells:list[str], hline:bool=False):
		assert len(cells)==self.n_colums, "number of colums is not correct! (" + str(self.n_colums) + " expected, " + str(len(cells)) + " given)"
	
		if self.text!="":
			self.text+= "\n"
			
		self.text+= " & ".join(cells)
		self.text+= r"\\"
		self.text+= "\n"
		
		if hline:
			self.text+= "\hline"
		

	def save(self):
		with open(self.destination, 'w') as target:
			target.write(self.text)

	def __str__(self):
		return "TexTable at " + self.destination + ":\n" + self.text
			
	 
