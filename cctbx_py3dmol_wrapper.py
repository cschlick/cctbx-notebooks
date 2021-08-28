from mmtbx.model.model import manager as model_manager
from iotbx.map_manager import map_manager
from iotbx.map_model_manager import map_model_manager
import py3Dmol
import numpy as np
import io


class CCTBX3dmolWrapper:
  """"
  A wrapper around the py3Dmol package written to 
  easily show CCTBX map/model objects.
  """

  def __init__(self,volume_cache_max=6):
    """
    3D volume are converted into a string before being send to 
    the 3dmol javascript library. This is slow, so we cache some 
    previously calculated volume strings.
    """
    self.volume_cache = {}
    self.volume_cache_stack = []
    self.volume_cache_max = volume_cache_max

  def clear_cache(self):
    self.volume_cache = {}
    self.volume_cache_stack = []

  def __call__(self,
              *objects,
              show_sidechains=True,
              show_ribbon=True,
              style=None,
              ribbon_color="red",
              sidechain_color=f"WhiteCarbon",
              map_color="#808080",
              opacity=0.6,
              map_threshold=0.1,
              map_apix=None,
              map_scale=1.0,
              map_shift=(0.0,0.0,0.0),
              debug=False):
    """
    The main function to show objects. Objects should be of type:
    
    mmtbx.model.model.manager
    iotbx.map_manager.map_manager
    iotbx.map_model_manager.map_model_manager

    When passing multiple objects at once, they will be displayed in a
    single viewport.

    Parameters:
    -----------
    objects: multiple cctbx objects to show
    show_sidechains (bool,True): Whether to show the protein sidechains
    show_ribbon (bool,True): Whether to show the protein ribbon
    style (dict,None): The py3Dmol style can be provided explicitly. 
                       ie: {'sphere':{"radius":1.0}} to display all atoms
                       This will overwrite the other display options.

    ribbon_color (str): The color of the ribbons
    sidechain_color (str): The color of the sidechains
    map_color (str): The color of the density map
    opacity (float): The opacity of the density map
    map_threshold (float): The density threshold
    map_apix (float): Angstroms per pixel. If None, taken from map_manager
    map_scale (float): modify map_apix by a scaler value
    map_shift (tuple): shift the map around



    """
    # collect all objects and sort as maps or models
    maps = []
    models = []
    for obj in objects:
      if isinstance(obj,map_manager):
        maps.append(obj)      
      elif isinstance(obj, map_model_manager):
        models.append(obj.model())
        maps.append(obj.map_manager())
      elif isinstance(obj,model_manager):
        models.append(obj)
      else:
        print("ERROR: Only high level CCTBX map/model objects are supported by this function")
        objs = [model_manager,map_manager,map_model_manager]
        for obj in objs:
          print("\n",obj.__class__)
        return None

    view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js',)

    ############# Maps #############
    map_strings = []
    for mm in maps:
      n_real = 1
      for n in mm.map_data().all():
        n_real*=n
      if n > 100**3:
        print("Map is large:",mm.map_data().all())
        print("This may take a very long time or crash")

      if map_apix is None:
        map_apix = mm.pixel_sizes()
      if hash(mm) in self.volume_cache.keys():
        cubestring = self.volume_cache[hash(mm)]
        if debug:
          print("Used volume cache:",True)
      else:
        data = mm.map_data().as_numpy_array()
        meta = self.mm_to_meta(mm)
        cubestring = write_cube_string(data,meta)
        if debug:
          print("Used volume cache:",False)
        self.volume_cache[hash(mm)] = cubestring
        self.volume_cache_stack.append(hash(mm))
        if len(self.volume_cache_stack)>self.volume_cache_max:
          del self.volume_cache[self.volume_cache_stack[0]]
          self.volume_cache_stack.pop(0)

      map_strings.append(cubestring)
    for map_string in map_strings:
      view.addVolumetricData(map_string, "cube", {'isoval': map_threshold, 'color':map_color , 'opacity': opacity})
    
    ############# Models #############
    for m in models:
      view.addModel(m.model_as_pdb(),'pdb')
    
    if style is not None:
       
       view.setStyle(style)
    else:
      if show_ribbon:
        view.setStyle({'cartoon': {'color':ribbon_color}})
      else:
        BB = ['C','O','N','CA']
        view.addStyle({'atom':BB},{'stick':{'colorscheme':sidechain_color,'radius':0.3}})
      if show_sidechains:
        BB = ['C','O','N']
        view.addStyle({'and':[{'resn':["GLY","PRO"],'invert':True},{'atom':BB,'invert':True}]},
                            {'stick':{'colorscheme':sidechain_color,'radius':0.3}})
        view.addStyle({'and':[{'resn':"GLY"},{'atom':'CA'}]},
                            {'sphere':{'colorscheme':sidechain_color,'radius':0.3}})
        view.addStyle({'and':[{'resn':"PRO"},{'atom':['C','O'],'invert':True}]},
                            {'stick':{'colorscheme':sidechain_color,'radius':0.3}})  
      
    view.zoomTo()
    return view
  

  @staticmethod
  def mm_to_meta(map_manager):
    """
    Form metadata dictionary from map manager object

    mm: map_manager objects

    returns: meta (dict)
    """
    conversion = 1.8897259885789233 # bohr to angstrom
    scale = conversion
    ax,ay,az = map_manager.pixel_sizes()
    ax*=scale
    ay*=scale
    az*=scale
    xvec = (ax,0.0,0.0)
    yvec = (0.0,ay,0.0)
    zvec = (0.0,0.0,az)
    origin = tuple(map_manager.get_origin())
    meta = {"atoms":[],
            "org":origin,
            "xvec":xvec,
            "yvec":yvec,
            "zvec":zvec}
    return meta
