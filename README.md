# FIRI-WHS

FIRI-WHS 农田水热盐数值模型
是一套模拟农田水分运动、热传导和盐分运移的命令行程序包，采用模块化程序设计，使用者可根据需求选用不同的模块进行模拟，也可修改现有模块或增加新的模块已满足不同模拟需求。

目前，我们完成了模型基本框架和部分水分运动模块的整理工作。
热传导、盐分运移、灌溉、作物、大气等模块的程序已基本实现，但是目前还需要时间进一步整理，工作完成将逐个发布。

土壤水分运动是围绕求解Richard方程构建的，需要输入土壤水分运动参数、初始条件、边界条件、源汇项等。

由于目前该程序包还在不断开发、修改和完善中，一些模块没有经过完全测试，请在使用过程中注意。

该程序包由中国农业科学院农田灌溉研究所黄仲冬、高阳，在洛桑研究所张效先的指导下进行开发。

联系方式: zdhuang@126.com

我们非常欢迎各位老师和同学批评指正，也欢迎进一步扩展该程序包。

在模型基本框架的搭建中，我们借鉴了MODFLOW早期版本的框架搭建思路，在程序的编写中，也采用了TTUTIL FORTRAN库，在次表示感谢！