# 图7至图11重绘说明

这五张图均使用R绘制，数据来自与当前手稿结果对应的已保存分析输出。没有重新拟合模型、筛选队列、修改患病率或调整R²。借鉴以下已发表论文的同类图形布局与视觉表达，再按本研究原有分面重新组合；并非复制论文图片或使用论文中的数据。组合图的每一个分面不一定都出现在同一参考图中。

| 本研究图号 | 主参考论文及具体图号 | 本次借鉴与保留内容 |
|---|---|---|
| 图7 | Lasky JR et al. Trait-mediated assembly processes predict successional changes in community diversity of tropical forests. PNAS. 2014;111:5616–5621. Figure 2. https://doi.org/10.1073/pnas.1319342111 | 横向排列标准化系数、简洁的点线标记。保留本研究的绝对系数及正负方向；没有借用参考图的区间。右侧保留本研究原有的顺序τ²展示。 |
| 图8 | Goldstein A et al. Peeking Inside the Black Box: Visualizing Statistical Learning With Plots of Individual Conditional Expectation. Journal of Computational and Graphical Statistics. 2015;24:44–65. Figure 3. https://doi.org/10.1080/10618600.2014.907095 | 淡色ICE曲线、突出显示平均预测、坐标轴上的数据分布刻线。保留本模型的真实阶梯形状和观察范围，不做平滑或老年外推。作者稿图3可核对：https://arxiv.org/pdf/1309.6392 。 |
| 图9 | Gehlenborg N, Wong B. Heat maps. Nature Methods. 2012;9:213. Figure 1. https://doi.org/10.1038/nmeth.1902 | 规则矩阵、统一色标和清楚的单元格边界。上下两部分仍分别为原有患病率热图与模型验证；负R²及重复验证范围原样保留。 |
| 图10 | A multi-ancestry polygenic risk score improves risk prediction for coronary artery disease. Nature Medicine. 2023. Figure 3, especially panels a, b and d. https://doi.org/10.1038/s41591-023-02429-x | 分层患病率点线图与分布箱线图的并排布局。仍使用本研究的Q1–Q4折外预测分层，未转为患者风险阈值。 |
| 图11 | Lever J, Krzywinski M, Altman N. Principal component analysis. Nature Methods. 2017;14:641–642. Figure 3. https://doi.org/10.1038/nmeth.4346 | PCA分组散点的简洁表达。保留原有雷达图，统一其配色、网格和文字；参考论文的Figure 3并不包含该雷达图。聚类坐标、组别和剖面数值均未重新估计。 |

通用样式为白底、Times New Roman、统一的加粗小写a/b/c标签，无冒号及阴影。蓝绿和朱红用于区分组别或数据类型。图7中增加样本量后的变化以“% of null”标示，相对于空模型τ²计算；它不是相对于前一阶段的相对增长率。图11中星号标记的患病率只用于聚类后的描述，没有进入聚类输入。

本仓库图件位于 `figures/SVG`，复现方法见根目录README。新版布局由 `code/Redraw_Figures_07_11.R` 生成；数据快照位于 `figure_data`。
