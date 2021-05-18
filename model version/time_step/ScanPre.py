#!/usr/bin/python3.6
# encoding:utf-8
from AddrsToSeq import InputAddrs
from Definitions import Stack
from DHC import SpaceTreeGen, OutputSpaceTree
from copy import deepcopy
import math
import pdb

def ScanPre(root,args):
    """
    动态扫描开始前的准备工作

    Args:
        root:空间树的根结点
    """

    InitializeDS(root, delta=args.delta)
    #InitializeTS(root) todo 不需要init了 每次需要多少就用多少就OK


def InitializeDS(node, parent_stack=Stack(), delta=16):
    """
    对结点node的DS进行初始化

    Args：
        node：当前DS待初始化的结点
        parent_stack：父结点的DS            
        delta：向量每一维度的基数
    """    
    
    # pdb.set_trace()
    parent=node.parent
    stack = deepcopy(parent_stack)  # 注意要将父结点的DS做拷贝
    if parent != None:
        stack.push(parent.diff_delta)
        node.DR.append(parent.diff_delta)
    vecDim = int(128 / math.log(delta, 2))
    for i in range(1, vecDim + 1):
        if node.Steady(i) and stack.find(i) == False:  # 如果这个位置same并且没有在父节点的stack中出现过就加入 szx
            stack.push(i)
            node.DR.append(i)
    if not node.isLeaf():
        for child in node.childs:
            InitializeDS(child, stack, delta)
    else:
        add_dim = {}
        for i in range(1, vecDim + 1):
            if stack.find(i) == False:
                add_dim[i] = node.get_entropy(i-1)
        add_dim = sorted(add_dim.items(), key=lambda kv: (kv[1], kv[0]))  # 实现按照熵值排序
        for i in add_dim:
            stack.push(i[0])
            node.DR.append(i[0])


    node.DS = stack

    node.R = len(node.iplist)/pow(delta, len(node.DR))     # 初始化 区域中种子数/未分配维度数 todo 有待修改


def InitializeTS(node):
    """
    对所有叶结点的TS进行初始化（SS和NDA在结点创建时已被初始化）

    Args：
        node：当前TS待初始化的结点
    """

    # pdb.set_trace()

    if node.isLeaf():
        delta = node.DS.pop()
        # print(node.node_id)
        # print(delta)
        # node.last_pop = delta
        # node.last_pop_value = node.TS[delta - 1]
        # print("leaf node :{}".format(node.global_node_id))
        node.ExpandTS(delta)
    else:
        for child in node.childs:
            InitializeTS(child)    
    # pdb.set_trace()
    


if __name__ == '__main__':
    IPS=InputAddrs(input="data1.csv")
    root=SpaceTreeGen(IPS,16,16)
    ScanPre(root)
    OutputSpaceTree(root)