#!/usr/bin/python3.6
# encoding:utf-8
import logging
import sys, os, platform
import nni
import argparse
from Definitions import Intersection, Stack, TreeNode
from AddrsToSeq import InputAddrs, SeqToAddrs
from DHC import SpaceTreeGen, OutputSpaceTree
from ScanPre import ScanPre
from ActiveScan import Scan
import pandas as pd
import numpy as np
from copy import deepcopy
import argparse
import math
import time
import random
if platform.system() == 'Linux':
    file = '/home/szx/code'
    sys.path.append(file)
    sys.path.extend([os.path.join(root, name) for root, dirs, _ in os.walk(file) for name in dirs])
#################################
"""
sudo python3 DynamicScan.py --input=/home/liguo/ipv6_project/6density/data1.csv --output=/home/liguo/ipv6_project/6density --budget=500  --IPv6=2001:da8:ff:212::10:3 --delta=16 --beta=16
"""
epochs = []
hit_rates = []
Remaining_times = []
scan_number = []
activate_number = []
def DynamicScan(root, args):
    """
    动态扫描空间树

    Args：
        root：空间树的根结点
        V：种子地址向量序列
        budget：扫描开销上限（最多扫描的地址数量）
        batch: 每轮扫描的次数
        alpha: 迭代参数
        source_ip：主机源IPv6地址
        output_dir：输出文件目录

    Return：
        R：经扫描发现的活跃地址集合（每个成员为真实的IPv6地址字符串）【序列？】
        P：检测到的别名前缀集合
        budget：剩余扫描次数
    """
    # OutputSpaceTree(root,V)
    # R = set()
    budget = args.budget
    batch = args.batch_size
    output_dir = args.output

    R = set()
    T = set()
    ALL_SAM = []
    now_budget = deepcopy(budget)
    active_file = output_dir + '/6density.result'+str(budget)
    target_file = output_dir + '/6density.target'+str(budget)
    xi = [] # 待扫描结点队列ξ
    InitializeNodeQueue(root, xi)
    last_budget = now_budget
    xi = sorted(xi, key=lambda node: node.R, reverse=True)
    epoch = 0
    print("节点数目：", len(xi))
    while 1:
        epochs.append(epoch)
        xi, now_budget, R, T = Scan_Feedback(xi, budget, now_budget, R, T, 0, args)
        if last_budget == now_budget:
            break
        else:
            last_budget = now_budget
        ReplaceDescendants(xi, args)
        xi = sorted(xi, key=lambda node: node.R, reverse=True)
        epoch+=1
        # 采用截取式 最终采样次数位于 [budget-batch, budget]
        if now_budget <= batch:
            #logging.info("Epoch: {}".format(epoch))
            xi, now_budget, R, T = Scan_Feedback(xi, budget, now_budget, R, T, 1, args)
            epochs.append(epoch)
            break
    dict_hit_result = {"epoch":epochs,"hit_rate":hit_rates,"remain_time":Remaining_times,"scaned_number":scan_number, "activate_number":activate_number}
    print(dict_hit_result)
    hit_result = pd.DataFrame(dict_hit_result)
    hit_result.to_csv(output_dir+'/budget_{},alpha={},num_node_p={},batch={},delta={},beta={}.csv'.format(args.budget, args.alpha, args.num_node_p, args.batch_size, args.delta, args.beta),index=False)
    print("begin to store the ip address in {} and {}".format(active_file, target_file))
    with open(active_file, 'w', encoding='utf-8') as f:
        for addr in R:
            f.write(addr + '\n')
    with open(target_file, 'w', encoding='utf-8') as f:
        for target in T:
            f.write(target + '\n')
    hit_rate = float(len(R))/(budget - now_budget)
    return R, budget-now_budget, len(R), hit_rate


def InitializeNodeQueue(root, xi):
    """
    层次遍历空间树，将结点队列ξ初始化为空间树的叶子结点

    Args：
        root:空间树的根结点
        xi：结点队列ξ
    """
    # pdb.set_trace()
    q = []
    q.append(root)
    while q != []:
        node = q.pop(0)
        if node.childs != []:  # 如果为空则为叶子节点
            q += node.childs
        else:
            xi.append(node)  # 存储叶子节点

def Scan_Feedback(xi, budget,now_budget, R, T,  end, args):
    """
    对队列xi中的所有结点进行一次扫描，
    并根据扫描得到的活跃地址密度对队列重新排序

    Args：
        xi：结点队列ξ
        budget：总扫描次数
        now_budget：剩余扫描次数
        batch_size：本次扫描的次数
        R：经扫描发现的活跃地址集合
        T
        V:种子地址向量集合
        source_ip
        output_dir
        target_file

    Return:
        xi:重新排序后的结点队列ξ
        budget:经过一次迭代扫描之后剩余的扫描次数
        R：更新后的活跃地址集合
        T：预测地址集合
    """

    # pdb.set_trace()

    batch_size = args.batch_size
    if now_budget < batch_size:
        batch_size = now_budget
    source_ip = args.IPv6
    output_dir = args.output
    alpha = args.alpha
    TS_addr_union = list()
    #temp_Sample_Set = []
    sum_Reward = 0.0
    max_reward = xi[0].R
    for i in range(max(int(len(xi)*args.num_node_p), 1)):  # 防止越界 todo 比例值防止等于0
        node = xi[i]
        sum_Reward += math.exp(node.R-max_reward)/(math.log(len(node.DR)+1))
    if sum_Reward<0.1:
        print("decrease node_num_p or increase batch_size")
        exit(1)

    for i in range(max(int(len(xi)*args.num_node_p), 1)):
        # if i % 100 == 0:
        #     print(i)
        node = xi[i]

        sample_result = SeqToAddrs(node, int(round((math.exp(node.R-max_reward)/math.log(len(node.DR)+1)/sum_Reward)*batch_size))+node.flag_child, args)
        #temp_Sample_Set.append((node, set(sample_result)))
        TS_addr_union += sample_result
        node.SS.update(set(sample_result))
        node.immediate_sample_size = len(sample_result)
    # for n in temp_Sample_Set:
    #     for m in ALL_SAM:
    #         if len(n[1].intersection(m[1])):
    #             print(n[1].intersection(m[1]))
    #             n[0].OutputNode()
    #             n[0].parent.OutputNode()
    #             m[0].OutputNode()
    #             m[0].parent.OutputNode()
    #             print("发现重复！")
    #             exit(1)
    # ALL_SAM += temp_Sample_Set
    C = set(TS_addr_union)
    now_budget -= len(C)
    T.update(C)
    scan_number.append(len(C))
    # with open(target_file, 'a', encoding='utf-8') as f:
    #     for target in C:
    #         f.write(target + '\n')
    active_addrs = set(Scan(C, source_ip, output_dir, 2))   #扫描并得到活跃的地址集合
    activate_number.append(len(active_addrs))
    print("intersection", len(R.intersection(active_addrs)))
    R.update(active_addrs)
    if end:
        nni.report_final_result(float(len(R) / (budget-now_budget)))
    else:
        nni.report_intermediate_result(float(len(R) / (budget-now_budget)))
    # 存储hit rate 与 remain_time 用于画图
    hit_rates.append(float(len(R) / (budget-now_budget)))
    Remaining_times.append(now_budget)
    print('[+]Hit rate:{}   Remaining scan times:{}\n'
       .format(float(len(R) / (budget-now_budget)), now_budget))

    for i in range(max(int(len(xi)*args.num_node_p), 1)):  # 更新节点的属性
        node = xi[i]
        new_active_addrs = active_addrs.intersection(node.SS)   # 寻找全部节点中与node.ss相似 的节点集合
        if node.immediate_sample_size > 0:
            node.R = (1-alpha)*node.R + alpha*float(len(new_active_addrs)/node.immediate_sample_size)  # 密度属性  更改成reward值
    xi = sorted(xi, key=lambda node: node.R, reverse=True)
    return xi, now_budget, R, T


def ReplaceDescendants(xi, args):
    """
    经过一次扫描后，若某结点的密度过大，在xi和xi_h队列中，
    需要将该结点及其所有的兄弟结点删除，并插入它们的父结点【证  明见Theorom3】

    Args：
        xi_h：下次将会被扫描的结点队列
        delta: 基数
    """
    # 将新节点加入集合中
    # 对新节点的子节点进行递归
    # 在递归的过程中合并 SS Reward  DR flag_child
    # 将子节点从中去除
    # pdb.set_trace()
    new_node = []
    for node in xi:
        if node.parent == None:
            continue   #注释！！！！！
        if node.parent.DS.stack == node.DS.stack:
            parent = node.parent
            print("begin merge {}".format(node.node_id))
            new_node.append(parent)  #用于最后判断是否已经删除了


            sum_node = Intersection(parent.all_childs, xi)

            sum_child_R = 0.0
            sum_all_DR = set(list(range(1,int(128 / math.log(args.delta, 2))+1)))

            sum_child_DR = sum_all_DR-set(parent.DS.stack)

            for node_i in sum_node:
                parent.SS.update(node_i.SS)
                sum_child_R += node_i.R*len(node_i.DR) # 合并R值
                parent.flag_child += node_i.flag_child
                # del node_i
            parent.searched_dim = len(sum_child_DR)
            parent.DR = list(set(parent.DR) | sum_child_DR)
            parent.R = sum_child_R/len(parent.DR) # 假定 父节点与子节点 空间大小与reward
            parent.ExpandTS(list(sum_child_DR))  #这一步的需要DR = list[]
            parent.region_size = pow(args.delta, parent.searched_dim) * len(parent.TS)
            for del_ in sum_node:
                try:
                    xi.remove(del_)  #防止其他节点未加入
                except:
                    pass
                try:
                    new_node.remove(del_)
                except:
                    pass

    for n in new_node:
        n.childs = []
        xi.append(n)


def MergeSort(xi_h, xi):
    """
    将两个有序的结点队列合并为一个

    Args：
        xi_h：队列1
        xi：队列2

    Return：
        queue：合并后的有序队列
    """

    queue = []
    i1 = 0
    i2 = 0
    while i1 < len(xi_h) or i2 < len(xi):
        if i1 >= len(xi_h):
            queue += xi[i2:]
            break
        elif i2 >= len(xi):
            queue += xi_h[i1:]
            break
        elif xi_h[i1].AAD >= xi[i2].AAD:
            queue.append(xi_h[i1])
            i1 += 1
        else:
            queue.append(xi[i2])
            i2 += 1

    return queue

def check(f, c, root):
    """
    用于判断两个节点是否有父子关系
    """
    q = []
    q.append(root)
    while q != []:
        node = q.pop(0)
        if node.node_id==c:
            fi = node
            break
        if node.childs != []:  # 如果为空则为叶子节点
            q += node.childs
    while fi.parent!=None:
        fi = fi.parent
        if fi.node_id == f:
            return True
    return False
def Start():
    parse=argparse.ArgumentParser()
    parse.add_argument('--input', type=str, default="/home/szx/code/baseline/V4/IPv6_dataset.txt",help='input IPv6 addresses')
    #parse.add_argument('--input', type=str, default="/home/szx/code/baseline/V4/testData.txt",help='input IPv6 addresses')
    parse.add_argument('--output',type=str, default="/home/szx/code/baseline/V4_alpha_down/result",help='output directory name')
    parse.add_argument('--budget',type=int, default=1000000, help='the upperbound of scan times')
    parse.add_argument('--IPv6',type=str, default="2001:da8:ff:212::7:8", help='local IPv6 address')
    parse.add_argument('--delta', type=int, default=4, help='the base of address')
    parse.add_argument('--beta',type=int, default=4, help='the max of node ')  # 1
    parse.add_argument('--alpha',type=float, default=1.0, help='the parameter of RL ')  # todo 设置成
    parse.add_argument('--batch_size', type=int, default=1000, help='the parameter for each epoch')  # todo 1. 每秒钟发现的活跃地址数 2. 命中率 设成维度倍数 小一点
    parse.add_argument('--num_node', type=int, default=10, help='the parameter for used nodes')     # 1 空间太大就不扫了
    parse.add_argument('--num_node_p', type=float, default=0.001, help='the parameter for used node probability')  # 1 空间太大就不扫了
    args_p=parse.parse_args()
    print(type(args_p))
    tuner_params = nni.get_next_parameter()
    args = vars(args_p)
    args.update(tuner_params)
    args = argparse.Namespace(**args)
    print(args)
    # args.input = '/home/sgl/6density_no_APD/files/source_copy.hex'
    # args.output = '/home/sgl/6density_no_APD/files2'
    # args.budget = 50000
    # args.IPv6 = '2001:da8:ff:212::20:5'

    # IPS=InputAddrs(input="data1.csv")
    # root=SpaceTreeGen(IPS,16,16)
    # OutputSpaceTree(root)
    #logging.basicConfig(filename='my_{},alpha={},num_node_p={},batch={},.log'.format(args.budget, args.alpha, args.num_node_p, args.batch_size), level=logging.INFO)
    print("ipv6 addres to sec begining")
    print(args.input)
    V = InputAddrs(input=args.input, delta=args.delta)
    print("ipv6 addres to sec over")
    print("SpaceTreeGen beginning")
    root = SpaceTreeGen(V, delta=args.delta, beta=args.beta)
    print("SpaceTreeGen over")
    ScanPre(root,args)
    # OutputSpaceTree(root,V)
    print('Space tree generated with {} seeds!'.format(len(V)))    
    R, target_len, result_len, hit_rate = DynamicScan(root, args)
    print('Over!')
    # hit_rate = float(len(R))/(init_budget - budget)
    # return init_budget - budget, len(R), hit_rate
    return target_len, result_len, hit_rate


if __name__ == '__main__':
    random.seed(1)
    target, result, hit_rate = Start()
    print("target {}".format(target))
    print("result {}".format(result))
    print("hit_rate {}".format(hit_rate))

    
