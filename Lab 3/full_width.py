"""Set jupyter notebook cellwidth to be wider"""
from IPython.core.display import display, HTML
def cell_width_percent(width_percent = '100'):
	"""Set jupyter notebook cellwidth to be wider"""
	assert (type(width_percent) == str), "width_percent must be string, 1->100"
	command = "<style>.container { width:" + width_percent +"% !important; }</style>"
	display(HTML(command))