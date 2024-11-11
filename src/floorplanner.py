###### SIMULATED ANNEALING ENGINE FOR SLICING FLOORPLAN ######

# Bash Terminal Commands Supported:

# [1] python3.7 python floorplanner.py -input HARD/n10.blocks -output n10.out  <-- Change source and output name as per .blocks file requirement
# Function: Provides a .out text file with processed coordinates of blocks with final and unused computed floorplanning area fetched from the SOFT/HARD .blocks File.

# Libraries used in this tool:
import argparse
import math
import random
from typing import List, Optional
import functools
HORIZONTAL_SLICE = -1
VERTICAL_SLICE = -2

# Assigned parameters for the SA engine:
STARTING_TEMP = 10000 # Setting the initial temperature
FREEZING_TEMP = 0.01 # Setting the final temperature
boltz_k = 1 # Boltzmann constant value

class WidthHeight:
    """
    Class declaration to map the object's Width and Height dimensionality.
    """
    def __init__(self, width: float = 0.0, height: float = 0.0, left_index: int = 0, right_index: int = 0):
        """
        Initializes the WidthHeight instance linked to parameters of width and height of the block
        """
        self.width = width # Value of Width
        self.height = height # Value of Height
        self.left_index = left_index # Index to represent Left boundary
        self.right_index = right_index # Index to represent Right boundary

    def rotate(self) -> 'WidthHeight':
        """
        Gives back a newly created instance from 'WidthHeight' with flipped dimensions.
        """
        return WidthHeight(self.height, self.width, self.left_index, self.right_index)

    def area(self) -> float:
        """
        Returns the area of the block by multiplying the width and height dimensions
        """
        return self.width * self.height

class Block:
    """
    Class declaration of a given block.
    """
    def __init__(self, block_name: str = "blank", isSoft: bool = False, widths_heights: Optional[List[WidthHeight]] = None,
                 leftChild: Optional['Block'] = None, rightChild: Optional['Block'] = None, parentBlock: Optional['Block'] = None):
        """
        Constructor method for Block. Checks whether the blocks contain Hard/Soft Macros.
        """
        self.block_name = block_name # Assigns name to the block
        self.isSoft = isSoft # Checker flag to verify whether the block is soft/hard macro.
        self.widths_heights = widths_heights or [] # Creates a list of WidthHeight objects representing dimensions
        self.leftChild = leftChild # Maps to the left child block
        self.rightChild = rightChild # Maps to the right child block
        self.parentBlock = parentBlock # Maps to the parent block
        self.wh_index = -1 # Index for widths_heights list
        self.xCoordinate = 0 # X-axis coordinate (Width) of the block
        self.yCoordinate = 0 # Y-axis coordinate (Height) of the block

    def rotate(self) -> 'Block':
        """
        Recreates again a new block with updated Width and Height dimensions.
        """
        new_wh = [wh.rotate() for wh in self.widths_heights]
        return Block(self.block_name, self.isSoft, new_wh)

    def deep_copy(self) -> 'Block':
        """
        Creates a deep copy of the block and its children.
        """
        left_child_copy = self.leftChild.deep_copy() if self.leftChild else None
        right_child_copy = self.rightChild.deep_copy() if self.rightChild else None
        copied_block = Block(self.block_name, self.isSoft,
                             [WidthHeight(wh.width, wh.height, wh.left_index, wh.right_index) for wh in self.widths_heights],
                             leftChild=left_child_copy, rightChild=right_child_copy)
        if left_child_copy:
            left_child_copy.parentBlock = copied_block
        if right_child_copy:
            right_child_copy.parentBlock = copied_block
        return copied_block
def parse_line_hard(line: str) -> Block:
    """
    Function to fetch hardrectilinear macro and assign it to the block object
    """
    parts = line.split()
    name = parts[0]
    width = float(parts[6].rstrip(',').strip(')'))
    height = float(parts[7].lstrip('(').rstrip(')').strip(','))
    wh = WidthHeight(width, height)
    rotated_wh = wh.rotate()
    return Block(name, isSoft=False, widths_heights=[wh, rotated_wh])

def parse_line_soft(line: str) -> Block:
    """
    Function to fetch softrectangular macro and assign it to the block object
    """
    parts = line.split()
    name = parts[0]
    area = float(parts[2])
    min_aspect = float(parts[3]) # Minimum aspect ratio of the block possible
    max_aspect = float(parts[4]) # Maximum aspect ratio of the block possible
    whs = []
    for aspect in (min_aspect, max_aspect):
        width = math.sqrt(area * aspect)
        height = math.sqrt(area / aspect)
        whs.append(WidthHeight(width, height))
        if min_aspect == max_aspect:
            break  # Avoids calculating duplicate blocks having aspect ratio.
    return Block(name, True, whs)


def parse_line(line):
    """
    Function to fetch line representing a hardrectilinear or softrectangular block and creates its respective Block object
    """
    split_line = line.split(' ')
    if split_line[1] == "hardrectilinear":
        return parse_line_hard(line)
    elif split_line[1] == "softrectangular":
        return parse_line_soft(line)


def read_blocks_from_file(filename):
    """
    Function for creating Block objects and reading block information from a file
    """
    blocks = []
    with open(filename, 'r') as fin:
        lines = fin.readlines()
        block_lines = lines[9:]
        for line in block_lines:

            if line.strip():
                block = parse_line(line)
                if block:
                    blocks.append(block)
    return blocks

class SlicingTree:
    """
    Class representing a slicing tree
    """
    def __init__(self, blks):
        self.blocks = []
        for b in blks:
            new_b = Block()  # Assuming Block has a duplicate or appropriate constructor
            new_b.__dict__ = b.__dict__.copy()
            self.blocks.append(new_b)
        for i in range(len(blks) - 1):
            if i % 2 == 1:
                self.blocks.append(Block("|"))
            else:
                self.blocks.append(Block("-"))

    def fix_pointers(self):
        """
        Function that fixes assigned pointers in the slicing tree
        """
        stack = []
        for i in self.blocks:
            if i.block_name == "-" or i.block_name == "|":
                assert len(stack) >= 2
                a = stack.pop()
                b = stack.pop()
                a.parent_block = i
                b.parent_block = i
                i.left_child = b
                i.right_child = a
            stack.append(i)

        assert len(stack) == 1
        stack[0].parentBlock = None
    def make_move(self):
        """
        Function to decide a move in the SA algorithm
        """
        self.score_up.cache_clear()
        self.score.cache_clear()
        valid_move = False
        new_score = None

        while not valid_move:
            move_type = random.choice(['orientation_change', 'swap','adjacent_swap'])

            if move_type == 'orientation_change':
                for rand_index, block in enumerate(self.blocks):
                    if block.block_name in "|-":
                        original_name = block.block_name
                        block.block_name = "-" if block.block_name == "|" else "|"
                        if self.is_valid():
                            new_score = self.score()
                            valid_move = True
                        else:
                            block.block_name = original_name
                        break

            elif move_type == 'swap':
                if len(self.blocks) >= 2:
                    r1, r2 = random.sample(range(len(self.blocks)), 2)
                    self.swap_blocks(r1, r2)
                    if self.is_valid():
                        new_score = self.score()
                        valid_move = True
                    else:
                        self.swap_blocks(r1, r2)  # Traceback

            elif move_type == 'adjacent_swap':
                candidates = [(i, i + 1) for i in range(len(self.blocks) - 1)
                              if (self.blocks[i].block_name in "|-" and self.blocks[i + 1].block_name not in "|-")
                              or (self.blocks[i].block_name not in "|-" and self.blocks[i + 1].block_name in "|-")]
                if candidates:
                    rand_index, adjacent_index = random.choice(candidates)
                    self.swap_blocks(rand_index, adjacent_index)
                    if self.is_valid():
                        new_score = self.score()
                        valid_move = True
                    else:
                        self.swap_blocks(rand_index, adjacent_index)  # Rollback

            self.fix_pointers()
            if valid_move:
                return new_score
        return None
    def swap_blocks(self, i, j):
        """
           Function to swap between 2 given blocks
        """
        
        self.blocks[i], self.blocks[j] = self.blocks[j], self.blocks[i]

    def is_valid(self):
        """
           Verifies if the slicing tree is valid.
        """
        count = 0
        prev_name = ""
        new_name = ""

        for b in self.blocks:
            new_name = b.block_name

            if new_name == "|" or new_name == "-":
                if count < 2:
                    return False
                if new_name == prev_name:
                    return False

                count -= 1
            else:
                count += 1

            prev_name = new_name

        return count == 1
    def score_single(self, b):
        """
           Computes the score of a single block.
        """
        if b.block_name == "|":
            assert b.left_child is not None
            assert b.right_child is not None

            vertical_node_sizing(b.left_child, b.right_child, b)
        elif b.block_name == "-":
            assert b.left_child is not None
            assert b.right_child is not None

            horizontal_node_sizing(b.left_child, b.right_child, b)
    @functools.lru_cache(maxsize=None)
    def score_up(self, b):
        """
           A Cached approach for computing a block's and its ancestors' scores.
        """
        self.score_single(b)
        if b.parentBlock is None:
            return min((wh.area() for wh in b.widths_heights), default=float('inf'))
        return self.score_up(b.parentBlock)

    @functools.lru_cache(maxsize=None)
    def score(self):
        """
           A Cached approach for computing the score of the slicing tree.
        """
        stack = []
        for i in self.blocks:
            if i.block_name == "-":
                assert len(stack) >= 2
                a = stack.pop()
                b = stack.pop()
                a.parent_block = i
                b.parent_block = i
                i.left_child = a
                i.right_child = b
                horizontal_node_sizing(a, b, i)
                stack.append(i)
            elif i.block_name == "|":
                assert len(stack) >= 2
                a = stack.pop()
                b = stack.pop()
                a.parent_block = i
                b.parent_block = i
                i.left_child = a
                i.right_child = b
                vertical_node_sizing(a, b, i)
                stack.append(i)
            else:
                stack.append(i)
        assert len(stack) == 1
        min_area = float('inf')
        min_area_index = -1

        assert len(stack[0].widths_heights) >= 1

        for i, wh in enumerate(stack[0].widths_heights):
            new_area = wh.area()

            if new_area < min_area:
                min_area = new_area
                min_area_index = i

        assert min_area_index != -1

        return min_area

    def copy_from(self, other):
        """
        Copies blocks from the other slicing tree.
        """
        self.blocks.clear()
        for block in other.blocks:
            self.blocks.append(block.deep_copy())

        self.fix_pointers()

    def compute_coords(self, b, wh_index):
        """
        Computes the coordinates of the blocks in the slicing tree.
        """
        if b.block_name == "|":
            assert b.left_child is not None and b.right_child is not None
            left_child_index = b.widths_heights[wh_index].left_index
            right_child_index = b.widths_heights[wh_index].right_index
            b.left_child.wh_index = left_child_index
            b.right_child.wh_index = right_child_index
            left_width = b.left_child.widths_heights[left_child_index].width
            b.left_child.xCoordinate = b.xCoordinate
            b.left_child.yCoordinate = b.yCoordinate
            b.right_child.xCoordinate = b.xCoordinate + left_width
            b.right_child.yCoordinate = b.yCoordinate
            self.compute_coords(b.left_child, left_child_index)
            self.compute_coords(b.right_child, right_child_index)

        elif b.block_name == "-":
            assert b.left_child is not None and b.right_child is not None
            left_child_index = b.widths_heights[wh_index].left_index
            right_child_index = b.widths_heights[wh_index].right_index
            b.left_child.wh_index = left_child_index
            b.right_child.wh_index = right_child_index
            left_height = b.left_child.widths_heights[left_child_index].height
            b.left_child.xCoordinate = b.xCoordinate
            b.left_child.yCoordinate = b.yCoordinate
            b.right_child.xCoordinate = b.xCoordinate
            b.right_child.yCoordinate = b.yCoordinate + left_height
            self.compute_coords(b.left_child, left_child_index)
            self.compute_coords(b.right_child, right_child_index)
        else:
            assert b.wh_index != -1

    def get_coords(self):
        """
        Fetches the coordinates from the blocks of the slicing tree as a string.
        """
        root = self.blocks[-1]

        assert root.block_name in ("|", "-")
        min_index = -1
        min_area = float('inf')

        assert len(root.widths_heights) >= 1

        for i, wh in enumerate(root.widths_heights):
            if wh.area() < min_area:
                min_area = wh.area()
                min_index = i

        assert min_index != -1
        self.compute_coords(root, min_index)
        coords_str = ""
        for b in self.blocks:
            if b.block_name not in ("|", "-"):
                assert b.wh_index != -1

                sx = b.xCoordinate
                sy = b.yCoordinate
                ex = sx + b.widths_heights[b.wh_index].width
                ey = sy + b.widths_heights[b.wh_index].height

                coords_str += f"{b.block_name} ({sx},{sy}) ({ex},{ey})\n"

        return coords_str

def floor_plan(blocks):
    """
    Performs floorplanning using the SA Engine.
    """
    random.seed(3)

    # Setting cooling factor based on count & type of the block.
    if any(block.isSoft for block in blocks) and len(blocks) == 300:
        NUM_MOVES_PER_STEP = 1.9437
        COOLING_FACTOR = 0.899
    elif any(block.isSoft for block in blocks) and len(blocks) == 200:
        NUM_MOVES_PER_STEP = 1.867
        COOLING_FACTOR = 0.899
    elif any(block.isSoft for block in blocks) and len(blocks) == 100:
        NUM_MOVES_PER_STEP = 1.55
        COOLING_FACTOR = 0.956
    elif any(not block.isSoft for block in blocks) and len(blocks) == 300:
        NUM_MOVES_PER_STEP = 1
        COOLING_FACTOR = 0.945
    else:
        NUM_MOVES_PER_STEP = 1
        COOLING_FACTOR = 0.927
        
        
    stree = SlicingTree(blocks)

    temp = STARTING_TEMP * len(blocks)
    freeze_temp = FREEZING_TEMP
    steps = int(NUM_MOVES_PER_STEP * len(blocks))
    current_score = stree.score()
    min_score = current_score
    min_stree = SlicingTree([])
    min_stree.copy_from(stree)
    min_stree.score_up.cache_clear()

    while temp > freeze_temp:
        for _ in range(steps):
            new_stree = SlicingTree([])
            new_stree.copy_from(stree)
            new_score = new_stree.make_move()
            assert new_stree.is_valid()
            cost = new_score - current_score
            if accept_move(cost, temp):
                stree.copy_from(new_stree)
                current_score = new_score
                if new_score < min_score:
                    min_score = new_score
                    min_stree.copy_from(new_stree)
        temp *= COOLING_FACTOR

    return min_stree

def accept_move(cost, temp):
    """
    Logic to decide whether to accept the move of block using SA algorithm
    """
    if cost < 0:
        return True
    else:
        boltz = math.exp(-1 * cost / (boltz_k * temp))
        r = random.random()  # Generates a random float number between 0.0 to 1.0
        return r < boltz

def sort_by_width(b):
    """
    Sorting blocks with respect to their individual widths.
    """
    b.widths_heights.sort(key=lambda wh: wh.width)

def sort_by_height(b):
    """
    Sorting blocks with respect to their individual heights.
    """
    b.widths_heights.sort(key=lambda wh: wh.height)
def vertical_node_sizing(a, b, result):
    """
    Calculation for sizing the vertical node.
    """
    sort_by_width(a)
    sort_by_width(b)
    result.widths_heights.clear()

    len_a = len(a.widths_heights)
    len_b = len(b.widths_heights)

    i, j = 0, 0
    while i < len_a and j < len_b:
        a_wh = a.widths_heights[i]
        b_wh = b.widths_heights[j]
        new_width = a_wh.width + b_wh.width
        new_height = max(a_wh.height, b_wh.height)
        new_wh = WidthHeight(width=new_width, height=new_height, left_index=i, right_index=j)
        result.widths_heights.append(new_wh)
        if new_height == a_wh.height:
            i += 1
        if new_height == b_wh.height:
            j += 1

def horizontal_node_sizing(a, b, result):
    """
    Calculation for sizing the horizontal node.
    """
    sort_by_height(a)
    sort_by_height(b)

    new_widths_heights = []  # Accumulate new WidthHeight instances here
    i, j = 0, 0
    len_a, len_b = len(a.widths_heights), len(b.widths_heights)

    while i < len_a and j < len_b:
        a_wh = a.widths_heights[i]
        b_wh = b.widths_heights[j]
        new_width = max(a_wh.width, b_wh.width)
        new_height = a_wh.height + b_wh.height
        new_wh = WidthHeight(width=new_width, height=new_height, left_index=i, right_index=j)
        new_widths_heights.append(new_wh)
        if new_width == a_wh.width:
            i += 1
        if new_width == b_wh.width:
            j += 1
    result.widths_heights = new_widths_heights

def output_solution(st, filename):
    """
    Outputs the processed solution to a file.
    """
    final_area = st.score()
    white_area = sum(b.widths_heights[0].area() for b in st.blocks if b.block_name not in ["|", "-"])
    with open(filename, 'w') as file:
        file.write(f"Final area = {final_area}\n")
        file.write(f"Black area = {final_area - white_area}\n")
        file.write("\n")
        file.write("block_name lower_left(x,y)coordinate upper_right(x,y)coordinate\n")
        file.write("\n")
        file.write(st.get_coords())  # Corrected method name

def main():
    """
    Main function and parsing command line arguments.
    """
    parser = argparse.ArgumentParser(description="Floorplanner")
    parser.add_argument("-input", help="Path to input file")
    parser.add_argument("-output", help="Path to output file")
    args = parser.parse_args()
    input_filename = args.input
    output_filename = args.output
    blocks = read_blocks_from_file(input_filename)
    optimized_tree= floor_plan(blocks)
    output_solution(optimized_tree, output_filename)


if __name__ == "__main__":
    main()
