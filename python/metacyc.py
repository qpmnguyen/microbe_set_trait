# Must have pathway tools open 
import pythoncyc

def get_hierarchy():
  meta = pythoncyc.select_organism("meta")
  path = meta.get_class_all_subs("|Pathways|")
  return(path)

def get_instances(class_name):
  meta = pythoncyc.select_organism("meta")
  inst = meta.get_class_data(class_name).instances
  inst_vec = [i.frameid for i in inst]
  return(inst_vec)
