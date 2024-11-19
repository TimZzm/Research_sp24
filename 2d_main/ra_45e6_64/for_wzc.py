# 把这个 DeDalus Plot 弄进来

import dedalus_Plot as DDP

# 第一步先创建一个Plot的instance，如果你不知道instance是什么，不用管，反正就是这么工作的
# 可以理解为你把你想要画的东西的数据存到了一个python variable里

# 你现在可以全弄默认值，就是括号里面不填东西
my_plot_1 = DDP.Plot()
print(f'\n\n\n\n\n')


# 也可以指定，注意save_dir可以直接相对路径而不用绝对路径。一般来说save_dir就是None，除非你真把数据存别的地方了
my_plot_2 = DDP.Plot(save_dir = None, handler = 'analysis_2')


"""
在这一步里，我们首先确定了我们的数据存放点，打开。
然后把所有的文件路径排了个序
再然后通过打开数据文件的最后一个文件来看你的h5文件的结构，并且打印出来方便查看
"""

# 下一步就可以画图了
# task name必填，这是你给你的物理变量的命名,
# 就是比如 analysis_2.add_task(LqW, name='liquid water')，在name= 那块的那个，你应该在上一步也能看到

# output_dir自己创建

# cmap看个人喜好，红蓝、viridis、spectral、rainbow都很实用，要是做水汽的话 Blues_r做出来的很贴近现实

# vmin, vmax如果有要求可以自己设，通常可以留空。
# 一般来说如果你要对比不同参数下某个物理变量的大小，可以手动设置它们，这样一来就可以直接用颜色对比两幅图里同一物理量大小。非常好的设计
# 这里我bobby做了一个改进就是保持了colorbar始终一致，用的是所有时刻的数据里平面出现的最大值和最小值，而不是每一时刻单独取平面最大最小值。
# 而且那个function的算法还能排查一下数据是否有大问题，比如出现NaN,inf的情况

# levelnum可以调高，就是你的colorbar分多少段，我喜欢32，图看起来基本等于连续了。

# figure_size可以完全随便的

# concentration 可以理解为广角镜头？我反正没怎么用过。我认为它的用处是可以调你想看的重点。我觉得 wzc现在就默认调1吧




my_plot_2.plot_all_snapshots(task_name='liquid water', output_dir='dir_for_zichu', cmap='viridis', 
                             vmin=None, vmax=None, levelnum=32, figure_size=(10, 8), concentration=1.0)


