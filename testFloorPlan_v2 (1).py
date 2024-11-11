##########################################################################################
# Permission to make digital or hard copies of all or part of this work for personal or
# classroom use is granted without fee provided that copies are not made or distributed
# for profit or commercial advantage and that copies bear this notice and the full citation
# on the first page. Copyrights for components of this work owned by others 
# must be honored. Abstracting with credit is permitted. To copy otherwise, or republish,
# to post on servers or to redistribute to lists, requires prior specific permission
# Request permissions from chhab011<at>umn<dot> edu. 
# For any inquires regarding this script , errors or change requests contact
# vachhabr<at>asu<dot>edu
#
#   Author:     Vidya Chhabria 
#   Email:      vachhabr@asu.edu
#   Revised:    15 - April - 2024 
#
#
##########################################################################################



import re
import sys

blocks = []

class Rect(object):
  def __init__(self, name, type, x2=None, y2=None, area=None, min=None, max=None):
    if type == 'hardrectilinear':
      self.name   = name
      self.type   = 'hardrectilinear'
      self.width  = x2
      self.length = y2
      self.area   = x2*y2
    elif type == 'softrectangular':
      self.name   = name
      self.type   = 'softrectangular'
      self.area   = area
      self.min    = min
      self.max    = max

  def update_location(self, x1, y1, x2, y2):
    self.left   = x1
    self.right  = x2
    self.bottom = y1
    self.top    = y2
    self.area1  = abs(x1-x2)*abs(y1-y2)

def read_input(file):
  print('Reading input file')
  error =0
  with open(file) as f:
    #for _ in range(8):
    #    next(f)
    #count = 9
    for line in f:
      if re.match('^[\t ]*[sS][bB]\d+',line):
        data = re.sub(r'[\(\),]','',line).split()
        if data[1] == 'hardrectilinear' and len(data) == 11:
          blocks.append(Rect(name = data[0], type='hardrectilinear', x2=float(data[7]), y2=float(data[8])))
        elif data[1] == 'softrectangular' and len(data) == 5:
          blocks.append(Rect(name=data[0], type='softrectangular', area=float(data[2]), min=float(data[3]), max=float(data[4])))
        else:
          print('ERROR - Input file at line number:', count, 'Unknown format')
          error =1
  if error == 0:
    print('CONFIRMED: Input file is in the correct format')
     #   count += 1
  
def read_output(file):
  def update_block(name, x1, y1, x2, y2):
    b = next((x for x in blocks if x.name == name), None)
    if b:
        b.update_location(x1, y1, x2, y2)
    else:
        print('\tERROR - block', name, 'is not found in input')

  reported_area_found = False
  black_area_found = False

  error =0
  print('Reading output file')
  #count = 4
  p = re.compile('^[+-]?[0-9]*\.?[0-9]*([Ee][+-]?[0-9]+)?$')
  pf = re.compile('[+-]?[\d]+\.?[\d]*([Ee][+-]?[0-9]+)?')
  with open(file) as f:
    for count, line in enumerate(f):
      if (not reported_area_found) and re.match('^[\t ]*[fF]inal [aA]rea[\t ]*=',line):
        #areas = re.findall('[+-]?[\d]+\.?[0-9]*([Ee][+-]?[0-9]+)?',line)
        m = pf.search(line)
        if m.group() is None:
          print(f'\tError - At line number {count+1} in output file. Only one area number was expected')
          error =1
        else:
          reported_area = abs(float(m.group()))
          reported_area_found =True
      if (not black_area_found) and re.match('^[\t ]*[bB]lack [aA]rea[\t ]*=',line):
        m = pf.search(line)
        if m.group() is None:
          print(f'\tError - At line number {count+1} in output file. Only one area number was expected')
          error =1
        else:
          black_area = abs(float(m.group()))
          black_area_found =True
       

      if re.match('^[\t ]*[sS][bB]\d+',line):
        data = re.findall(r'[\w\.\-\+]+', line)
        coord = [x for x in data if p.match(x)]
        if len(data) != 5:
          print('\tError - At line number', count+1, 'in output file. 5 words expected')
          error =1
        elif len(coord) != 4:
          print('\tError - At line number', count+1, 'in output file. 4 numbers expected')
          error =1
        else:
          update_block(data[0], float(data[1]), float(data[2]), float(data[3]), float(data[4]))
  if (not reported_area_found):
    print("Error - The final area of the block is not reported")
  if (not black_area_found):
    print("Error - The black area of the block is not reported")
  if error == 0:
    print('CONFIRMED: Output file is in the correct format')
  return reported_area, black_area

def check_area(reported_area, black_area):
  error =False
  max_x, max_y = float("-inf"), float("-inf") 
  min_x, min_y = float("inf"), float("inf") 
  cumulative_area = 0
  for b in blocks:
    min_x = min( min_x, b.left, b.right )
    max_x = max( max_x, b.left, b.right )
    min_y = min( min_y, b.top,  b.bottom)
    max_y = max( max_y, b.top,  b.bottom)
    #print(b.name, b.left, b.right, b.top, b.bottom)
    cumulative_area += (abs(b.left - b.right)*abs(b.top-b.bottom))
    
  calculated_area = (max_x-min_x)*(max_y-min_y)
  calculated_black_area = calculated_area - cumulative_area
  if (abs(calculated_area - reported_area)>0.01*reported_area):
    print(f"Error - Reported final area does not match the floor plan")
    print(f"Calculated Area: {calculated_area}, Reported: {reported_area}")
    print(f"overall coordinates: x1: {min_x}, x2: {max_x}, y1:{min_y}, y2: {max_y}")
    error=True
  if(abs(calculated_black_area - black_area)>0.01*black_area):
    print(f"Error - Reported black area does not match the floor plan")
    print(f"Cumulative black area computed: {calculated_black_area}, reported: {black_area}")
    error=True
  if(calculated_black_area > 0.25*calculated_area):
    print(f"Error - >25% empty area in the chip")
    error=True
  if not error:
    print('CONFIRMED: All areas are validated')


def check_missing():
    print('Checking area mismatch and missing blocks')
    error =0
    for b in blocks:
        if not hasattr(b, 'area1'):
            print('\tERROR - block', b.name, 'is not found in output')
            error =1
        elif abs(b.area1 - b.area) > 0.01*b.area:
            print('\tERROR - For block', b.name, 'area is not matching')
            error =1
    if error == 0:
        print('CONFIRMED: All blocks accounted for with correct areas ')


def check_overlaps():
    def overlap(r1, r2):
        h_overlaps = (r1.left < r2.right) and (r1.right > r2.left)
        v_overlaps = (r1.bottom < r2.top) and (r1.top > r2.bottom)
        return h_overlaps and v_overlaps

    error =0
    print('Checking overlapping blocks')
    for b1 in blocks:
        for b2 in blocks:
            if b1.name == b2.name:
                continue
            elif overlap(b1, b2):
                print('\tERROR -', b1.name, 'and', b2.name, 'overlap')
                print('\nThere might be more overlapping blocks, but displaying only one here')
                print('Please correct it and test again')
                error =1
                sys.exit()
    if error == 0:
        print('CONFIRMED: No overlapping blocks detected')

def check_aspect_ratio():
  error =0
  print('Checking rotation for hard-macros')
  print('Checking aspect ratio for soft-macros')
  for b in blocks:
    if b.type == 'softrectangular':
      width = abs(b.right - b.left)
      length = abs(b.top - b.bottom)
      difference_max = abs(b.max - width/length)
      difference_min = abs(b.min - width/length)
      difference_same = abs(1 - width/length)
      if (1<=b.max) and (1>=b.min):
        if (difference_max > 0.01) and (difference_min > 0.01) and (difference_same > 0.001):
          print('\tERROR -', b.name, 'has wrong aspect ratio')
          error =1
      else:
        if (difference_max > 0.01) and (difference_min > 0.01):
          print('\tERROR -', b.name, 'has wrong aspect ratio')
          error =1
    else: 
      reported_width = abs(b.right - b.left)
      reported_length = abs(b.top - b.bottom)
      if (not (((abs(reported_length - b.length)<0.01*b.length)
            and (abs(reported_width - b.width)<0.01*b.width))
            or ((abs(reported_length - b.width)<0.01*b.width)
            and (abs(reported_width - b.length)<0.01*b.length))
         )):
        print(f'\tError - {b.name} is does not have the correct rotation')
        

  if error == 0:
      print('CONFIRMED: Aspect ratios are in-order')

def display_header():
  print('\n\
  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n\
  +  _____   _        ___     ___    ____    ____    _          _      _   _        _____   _____   ____    _____   _____   ____    + \n\
  + |  ___| | |      / _ \   / _ \  |  _ \  |  _ \  | |        / \    | \ | |      |_   _| | ____| / ___|  |_   _| | ____| |  _ \   + \n\
  + | |_    | |     | | | | | | | | | |_) | | |_) | | |       / _ \   |  \| |        | |   |  _|   \___ \    | |   |  _|   | |_) |  + \n\
  + |  _|   | |___  | |_| | | |_| | |  _ <  |  __/  | |___   / ___ \  | |\  |        | |   | |___   ___) |   | |   | |___  |  _ <   + \n\
  + |_|     |_____|  \___/   \___/  |_| \_\ |_|     |_____| /_/   \_\ |_| \_|        |_|   |_____| |____/    |_|   |_____| |_| \_\  + \n\
  +                                                                                                                                 + \n\
  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n\
  ')



def main():
  if len(sys.argv) == 3:
    display_header();
    read_input(sys.argv[1])
    reported_area, black_area = read_output(sys.argv[2])
    check_missing()
    check_area(reported_area, black_area)
    check_overlaps()
    check_aspect_ratio()
  else:
    print('USAGE: ./testFloorPlan.py inputFileName outputFileName')

main()