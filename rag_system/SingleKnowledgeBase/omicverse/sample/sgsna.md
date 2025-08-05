# sgsna
*Converted from: omicverse/sample/sgsna.ipynb*

```pythonimport networkx as nx
G_all=nx.Graph()
for i in c:
  G_all.add_edge(i[0],i[1])```

```pythonscipy_kde=st.gaussian_kde(data)
dens = scipy_kde(data)
dens```
*Output:*

```pythonimport numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
plt.figure(figsize = (10, 8))
X = np.array(data['F2_3'])  # 转化为1D array
scipy_kde = st.gaussian_kde(X)  # 高斯核密度估计
X.sort()
'''这三种方法都可以得到估计的概率密度'''
dens = scipy_kde.evaluate(X)
#dens1=scipy_kde.pdf(X)  # pdf求概率密度
#dens2=scipy_kde(X)
plt.plot(X, dens, c='green', label='核密度值')
plt.tick_params(labelsize = 20)  # 设置坐标刻度值的大小      
font = {'size': 20}  # 设置横纵坐标的名称以及对应字体格式、大小
plt.xlabel('变量', font)
plt.ylabel('概率密度函数', font)
plt.legend(fontsize = 15)  # 显示图例,设置图例字体大小
plt.show()```
*Output:*
```<Figure size 720x576 with 1 Axes>```

```pythondens```
*Output:*
```array([0.00497294, 0.0115484 , 0.01253968, ..., 0.02311075, 0.01758854,
       0.01138771])```

# 软阈值β筛选

## 1、导入数据

```pythonimport pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
from scipy.stats import norm
from scipy import stats
import networkx as nx
import datetime

data=pd.read_csv('LiverFemale3600.csv')
data.dropna(inplace=True)
data.set_index(data.columns[0],inplace=True)
data.head()```
*Output:*
```                F2_2    F2_3     F2_14    F2_15    F2_19     F2_20    F2_23  \
substanceBXH                                                                  
MMT00000044  -0.0181  0.0642  0.000064 -0.05800  0.04830 -0.151974 -0.00129   
MMT00000046  -0.0773 -0.0297  0.112000 -0.05890  0.04430 -0.093800  0.09340   
MMT00000051  -0.0226  0.0617 -0.129000  0.08710 -0.11500 -0.065026  0.00249   
MMT00000080  -0.0487  0.0582 -0.048300 -0.03710  0.02510  0.085043  0.04450   
MMT00000102   0.1760 -0.1890 -0.065000 -0.00846 -0.00574 -0.018072 -0.12500   

                F2_24   F2_26    F2_37    ...       F2_324  F2_325  F2_326  \
substanceBXH                              ...                                
MMT00000044  -0.23600 -0.0307 -0.02610    ...     0.047700 -0.0488  0.0168   
MMT00000046   0.02690 -0.1330  0.07570    ...    -0.049200 -0.0350 -0.0738   
MMT00000051  -0.10200  0.1420 -0.10200    ...     0.000612  0.1210  0.0996   
MMT00000080   0.00167 -0.0680  0.00567    ...     0.113000 -0.0859 -0.1340   
MMT00000102  -0.06820  0.1250  0.00998    ...    -0.080000 -0.1200  0.1230   

              F2_327   F2_328  F2_329  F2_330  F2_332  F2_355    F2_357  
substanceBXH                                                             
MMT00000044  -0.0309  0.02740 -0.0310  0.0660 -0.0199 -0.0146  0.065000  
MMT00000046  -0.1730 -0.07380 -0.2010 -0.0820 -0.0939  0.0192 -0.049900  
MMT00000051   0.1090  0.02730  0.1200 -0.0629 -0.0395  0.1090  0.000253  
MMT00000080   0.0639  0.00731  0.1240 -0.0212  0.0870  0.0512  0.024300  
MMT00000102   0.1870  0.05410  0.0699  0.0708  0.1450 -0.0399  0.037500  

[5 rows x 135 columns]```

```pythonhelp(data.T.corr())```
*Output:*
```Help on DataFrame in module pandas.core.frame object:

class DataFrame(pandas.core.generic.NDFrame)
 |  DataFrame(data=None, index=None, columns=None, dtype=None, copy=False)
 |  
 |  Two-dimensional size-mutable, potentially heterogeneous tabular data
 |  structure with labeled axes (rows and columns). Arithmetic operations
 |  align on both row and column labels. Can be thought of as a dict-like
 |  container for Series objects. The primary pandas data structure.
 |  
 |  Parameters
 |  ----------
 |  data : numpy ndarray (structured or homogeneous), dict, or DataFrame
 |      Dict can contain Series, arrays, constants, or list-like objects
 |  
 |      .. versionchanged :: 0.23.0
 |         If data is a dict, argument order is maintained for Python 3.6
 |         and later.
 |  
 |  index : Index or array-like
 |      Index to use for resulting frame. Will default to RangeIndex if
 |      no indexing information part of input data and no index provided
 |  columns : Index or array-like
 |      Column labels to use for resulting frame. Will default to
 |      RangeIndex (0, 1, 2, ..., n) if no column labels are provided
 |  dtype : dtype, default None
 |      Data type to force. Only a single dtype is allowed. If None, infer
 |  copy : boolean, default False
 |      Copy data from inputs. Only affects DataFrame / 2d ndarray input
 |  
 |  Examples
 |  --------
 |  Constructing DataFrame from a dictionary.
 |  
 |  >>> d = {'col1': [1, 2], 'col2': [3, 4]}
 |  >>> df = pd.DataFrame(data=d)
 |  >>> df
 |     col1  col2
 |  0     1     3
 |  1     2     4
 |  
 |  Notice that the inferred dtype is int64.
 |  
 |  >>> df.dtypes
 |  col1    int64
 |  col2    int64
 |  dtype: object
 |  
 |  To enforce a single dtype:
 |  
 |  >>> df = pd.DataFrame(data=d, dtype=np.int8)
 |  >>> df.dtypes
 |  col1    int8
 |  col2    int8
 |  dtype: object
 |  
 |  Constructing DataFrame from numpy ndarray:
 |  
 |  >>> df2 = pd.DataFrame(np.random.randint(low=0, high=10, size=(5, 5)),
 |  ...                    columns=['a', 'b', 'c', 'd', 'e'])
 |  >>> df2
 |      a   b   c   d   e
 |  0   2   8   8   3   4
 |  1   4   2   9   0   9
 |  2   1   0   7   8   0
 |  3   5   1   7   1   3
 |  4   6   0   2   4   2
 |  
 |  See also
 |  --------
 |  DataFrame.from_records : constructor from tuples, also record arrays
 |  DataFrame.from_dict : from dicts of Series, arrays, or dicts
 |  DataFrame.from_items : from sequence of (key, value) pairs
 |  pandas.read_csv, pandas.read_table, pandas.read_clipboard
 |  
 |  Method resolution order:
 |      DataFrame
 |      pandas.core.generic.NDFrame
 |      pandas.core.base.PandasObject
 |      pandas.core.base.StringMixin
 |      pandas.core.accessor.DirNamesMixin
 |      pandas.core.base.SelectionMixin
 |      builtins.object
 |  
 |  Methods defined here:
 |  
 |  __add__(self, other, axis=None, level=None, fill_value=None)
 |      Binary operator __add__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __and__(self, other, axis='columns', level=None, fill_value=None)
 |      Binary operator __and__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __div__ = __truediv__(self, other, axis=None, level=None, fill_value=None)
 |  
 |  __eq__(self, other)
 |      Wrapper for comparison method __eq__
 |  
 |  __floordiv__(self, other, axis=None, level=None, fill_value=None)
 |      Binary operator __floordiv__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __ge__(self, other)
 |      Wrapper for comparison method __ge__
 |  
 |  __getitem__(self, key)
 |  
 |  __gt__(self, other)
 |      Wrapper for comparison method __gt__
 |  
 |  __iadd__ = f(self, other)
 |  
 |  __iand__ = f(self, other)
 |  
 |  __ifloordiv__ = f(self, other)
 |  
 |  __imod__ = f(self, other)
 |  
 |  __imul__ = f(self, other)
 |  
 |  __init__(self, data=None, index=None, columns=None, dtype=None, copy=False)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  __ior__ = f(self, other)
 |  
 |  __ipow__ = f(self, other)
 |  
 |  __isub__ = f(self, other)
 |  
 |  __itruediv__ = f(self, other)
 |  
 |  __ixor__ = f(self, other)
 |  
 |  __le__(self, other)
 |      Wrapper for comparison method __le__
 |  
 |  __len__(self)
 |      Returns length of info axis, but here we use the index
 |  
 |  __lt__(self, other)
 |      Wrapper for comparison method __lt__
 |  
 |  __matmul__(self, other)
 |      Matrix multiplication using binary `@` operator in Python>=3.5
 |  
 |  __mod__(self, other, axis=None, level=None, fill_value=None)
 |      Binary operator __mod__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __mul__(self, other, axis=None, level=None, fill_value=None)
 |      Binary operator __mul__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __ne__(self, other)
 |      Wrapper for comparison method __ne__
 |  
 |  __or__(self, other, axis='columns', level=None, fill_value=None)
 |      Binary operator __or__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __pow__(self, other, axis=None, level=None, fill_value=None)
 |      Binary operator __pow__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __radd__(self, other, axis=None, level=None, fill_value=None)
 |      Binary operator __radd__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __rand__(self, other, axis='columns', level=None, fill_value=None)
 |      Binary operator __rand__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __rdiv__ = __rtruediv__(self, other, axis=None, level=None, fill_value=None)
 |  
 |  __rfloordiv__(self, other, axis=None, level=None, fill_value=None)
 |      Binary operator __rfloordiv__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __rmatmul__(self, other)
 |      Matrix multiplication using binary `@` operator in Python>=3.5
 |  
 |  __rmod__(self, other, axis=None, level=None, fill_value=None)
 |      Binary operator __rmod__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __rmul__(self, other, axis=None, level=None, fill_value=None)
 |      Binary operator __rmul__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __ror__(self, other, axis='columns', level=None, fill_value=None)
 |      Binary operator __ror__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __rpow__(self, other, axis=None, level=None, fill_value=None)
 |      Binary operator __rpow__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __rsub__(self, other, axis=None, level=None, fill_value=None)
 |      Binary operator __rsub__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __rtruediv__(self, other, axis=None, level=None, fill_value=None)
 |      Binary operator __rtruediv__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __rxor__(self, other, axis='columns', level=None, fill_value=None)
 |      Binary operator __rxor__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __setitem__(self, key, value)
 |  
 |  __sub__(self, other, axis=None, level=None, fill_value=None)
 |      Binary operator __sub__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __truediv__(self, other, axis=None, level=None, fill_value=None)
 |      Binary operator __truediv__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  __unicode__(self)
 |      Return a string representation for a particular DataFrame
 |      
 |      Invoked by unicode(df) in py2 only. Yields a Unicode String in both
 |      py2/py3.
 |  
 |  __xor__(self, other, axis='columns', level=None, fill_value=None)
 |      Binary operator __xor__ with support to substitute a fill_value for missing data in
 |      one of the inputs
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |  
 |  add(self, other, axis='columns', level=None, fill_value=None)
 |      Addition of dataframe and other, element-wise (binary operator `add`).
 |      
 |      Equivalent to ``dataframe + other``, but with support to substitute a fill_value for
 |      missing data in one of the inputs.
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      
 |      >>> a = pd.DataFrame([1, 1, 1, np.nan], index=['a', 'b', 'c', 'd'],
 |      ...                  columns=['one'])
 |      >>> a
 |         one
 |      a  1.0
 |      b  1.0
 |      c  1.0
 |      d  NaN
 |      >>> b = pd.DataFrame(dict(one=[1, np.nan, 1, np.nan],
 |      ...                       two=[np.nan, 2, np.nan, 2]),
 |      ...                  index=['a', 'b', 'd', 'e'])
 |      >>> b
 |         one  two
 |      a  1.0  NaN
 |      b  NaN  2.0
 |      d  1.0  NaN
 |      e  NaN  2.0
 |      >>> a.add(b, fill_value=0)
 |         one  two
 |      a  2.0  NaN
 |      b  1.0  2.0
 |      c  1.0  NaN
 |      d  1.0  NaN
 |      e  NaN  2.0
 |      
 |      
 |      See also
 |      --------
 |      DataFrame.radd
 |  
 |  agg = aggregate(self, func, axis=0, *args, **kwargs)
 |  
 |  aggregate(self, func, axis=0, *args, **kwargs)
 |      Aggregate using one or more operations over the specified axis.
 |      
 |      .. versionadded:: 0.20.0
 |      
 |      Parameters
 |      ----------
 |      func : function, string, dictionary, or list of string/functions
 |          Function to use for aggregating the data. If a function, must either
 |          work when passed a DataFrame or when passed to DataFrame.apply. For
 |          a DataFrame, can pass a dict, if the keys are DataFrame column names.
 |      
 |          Accepted combinations are:
 |      
 |          - string function name.
 |          - function.
 |          - list of functions.
 |          - dict of column names -> functions (or list of functions).
 |      
 |      
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          - 0 or 'index': apply function to each column.
 |          - 1 or 'columns': apply function to each row.
 |      *args
 |          Positional arguments to pass to `func`.
 |      **kwargs
 |          Keyword arguments to pass to `func`.
 |      
 |      Returns
 |      -------
 |      aggregated : DataFrame
 |      
 |      Notes
 |      -----
 |      `agg` is an alias for `aggregate`. Use the alias.
 |      
 |      A passed user-defined-function will be passed a Series for evaluation.
 |      
 |      The aggregation operations are always performed over an axis, either the
 |      index (default) or the column axis. This behavior is different from
 |      `numpy` aggregation functions (`mean`, `median`, `prod`, `sum`, `std`,
 |      `var`), where the default is to compute the aggregation of the flattened
 |      array, e.g., ``numpy.mean(arr_2d)`` as opposed to ``numpy.mean(arr_2d,
 |      axis=0)``.
 |      
 |      `agg` is an alias for `aggregate`. Use the alias.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame([[1, 2, 3],
 |      ...                    [4, 5, 6],
 |      ...                    [7, 8, 9],
 |      ...                    [np.nan, np.nan, np.nan]],
 |      ...                   columns=['A', 'B', 'C'])
 |      
 |      Aggregate these functions over the rows.
 |      
 |      >>> df.agg(['sum', 'min'])
 |              A     B     C
 |      sum  12.0  15.0  18.0
 |      min   1.0   2.0   3.0
 |      
 |      Different aggregations per column.
 |      
 |      >>> df.agg({'A' : ['sum', 'min'], 'B' : ['min', 'max']})
 |              A    B
 |      max   NaN  8.0
 |      min   1.0  2.0
 |      sum  12.0  NaN
 |      
 |      Aggregate over the columns.
 |      
 |      >>> df.agg("mean", axis="columns")
 |      0    2.0
 |      1    5.0
 |      2    8.0
 |      3    NaN
 |      dtype: float64
 |      
 |      See also
 |      --------
 |      DataFrame.apply : Perform any type of operations.
 |      DataFrame.transform : Perform transformation type operations.
 |      pandas.core.groupby.GroupBy : Perform operations over groups.
 |      pandas.core.resample.Resampler : Perform operations over resampled bins.
 |      pandas.core.window.Rolling : Perform operations over rolling window.
 |      pandas.core.window.Expanding : Perform operations over expanding window.
 |      pandas.core.window.EWM : Perform operation over exponential weighted
 |          window.
 |  
 |  align(self, other, join='outer', axis=None, level=None, copy=True, fill_value=None, method=None, limit=None, fill_axis=0, broadcast_axis=None)
 |      Align two objects on their axes with the
 |      specified join method for each axis Index
 |      
 |      Parameters
 |      ----------
 |      other : DataFrame or Series
 |      join : {'outer', 'inner', 'left', 'right'}, default 'outer'
 |      axis : allowed axis of the other object, default None
 |          Align on index (0), columns (1), or both (None)
 |      level : int or level name, default None
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      copy : boolean, default True
 |          Always returns new objects. If copy=False and no reindexing is
 |          required then original objects are returned.
 |      fill_value : scalar, default np.NaN
 |          Value to use for missing values. Defaults to NaN, but can be any
 |          "compatible" value
 |      method : str, default None
 |      limit : int, default None
 |      fill_axis : {0 or 'index', 1 or 'columns'}, default 0
 |          Filling axis, method and limit
 |      broadcast_axis : {0 or 'index', 1 or 'columns'}, default None
 |          Broadcast values along this axis, if aligning two objects of
 |          different dimensions
 |      
 |      Returns
 |      -------
 |      (left, right) : (DataFrame, type of other)
 |          Aligned objects
 |  
 |  all(self, axis=0, bool_only=None, skipna=True, level=None, **kwargs)
 |      Return whether all elements are True, potentially over an axis.
 |      
 |      Returns True if all elements within a series or along a Dataframe
 |      axis are non-zero, not-empty or not-False.
 |      
 |      Parameters
 |      ----------
 |      axis : {0 or 'index', 1 or 'columns', None}, default 0
 |          Indicate which axis or axes should be reduced.
 |      
 |          * 0 / 'index' : reduce the index, return a Series whose index is the
 |            original column labels.
 |          * 1 / 'columns' : reduce the columns, return a Series whose index is the
 |            original index.
 |          * None : reduce all axes, return a scalar.
 |      
 |      skipna : boolean, default True
 |          Exclude NA/null values. If an entire row/column is NA, the result
 |          will be NA.
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series.
 |      bool_only : boolean, default None
 |          Include only boolean columns. If None, will attempt to use everything,
 |          then use only boolean data. Not implemented for Series.
 |      **kwargs : any, default None
 |          Additional keywords have no effect but might be accepted for
 |          compatibility with NumPy.
 |      
 |      Returns
 |      -------
 |      all : Series or DataFrame (if level specified)
 |      
 |      See also
 |      --------
 |      pandas.Series.all : Return True if all elements are True
 |      pandas.DataFrame.any : Return True if one (or more) elements are True
 |      
 |      Examples
 |      --------
 |      Series
 |      
 |      >>> pd.Series([True, True]).all()
 |      True
 |      >>> pd.Series([True, False]).all()
 |      False
 |      
 |      DataFrames
 |      
 |      Create a dataframe from a dictionary.
 |      
 |      >>> df = pd.DataFrame({'col1': [True, True], 'col2': [True, False]})
 |      >>> df
 |         col1   col2
 |      0  True   True
 |      1  True  False
 |      
 |      Default behaviour checks if column-wise values all return True.
 |      
 |      >>> df.all()
 |      col1     True
 |      col2    False
 |      dtype: bool
 |      
 |      Specify ``axis='columns'`` to check if row-wise values all return True.
 |      
 |      >>> df.all(axis='columns')
 |      0     True
 |      1    False
 |      dtype: bool
 |      
 |      Or ``axis=None`` for whether every value is True.
 |      
 |      >>> df.all(axis=None)
 |      False
 |  
 |  any(self, axis=0, bool_only=None, skipna=True, level=None, **kwargs)
 |      Return whether any element is True over requested axis.
 |      
 |      Unlike :meth:`DataFrame.all`, this performs an *or* operation. If any of the
 |      values along the specified axis is True, this will return True.
 |      
 |      Parameters
 |      ----------
 |      axis : {0 or 'index', 1 or 'columns', None}, default 0
 |          Indicate which axis or axes should be reduced.
 |      
 |          * 0 / 'index' : reduce the index, return a Series whose index is the
 |            original column labels.
 |          * 1 / 'columns' : reduce the columns, return a Series whose index is the
 |            original index.
 |          * None : reduce all axes, return a scalar.
 |      
 |      skipna : boolean, default True
 |          Exclude NA/null values. If an entire row/column is NA, the result
 |          will be NA.
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series.
 |      bool_only : boolean, default None
 |          Include only boolean columns. If None, will attempt to use everything,
 |          then use only boolean data. Not implemented for Series.
 |      **kwargs : any, default None
 |          Additional keywords have no effect but might be accepted for
 |          compatibility with NumPy.
 |      
 |      Returns
 |      -------
 |      any : Series or DataFrame (if level specified)
 |      
 |      See Also
 |      --------
 |      pandas.DataFrame.all : Return whether all elements are True.
 |      
 |      Examples
 |      --------
 |      **Series**
 |      
 |      For Series input, the output is a scalar indicating whether any element
 |      is True.
 |      
 |      >>> pd.Series([True, False]).any()
 |      True
 |      
 |      **DataFrame**
 |      
 |      Whether each column contains at least one True element (the default).
 |      
 |      >>> df = pd.DataFrame({"A": [1, 2], "B": [0, 2], "C": [0, 0]})
 |      >>> df
 |         A  B  C
 |      0  1  0  0
 |      1  2  2  0
 |      
 |      >>> df.any()
 |      A     True
 |      B     True
 |      C    False
 |      dtype: bool
 |      
 |      Aggregating over the columns.
 |      
 |      >>> df = pd.DataFrame({"A": [True, False], "B": [1, 2]})
 |      >>> df
 |             A  B
 |      0   True  1
 |      1  False  2
 |      
 |      >>> df.any(axis='columns')
 |      0    True
 |      1    True
 |      dtype: bool
 |      
 |      >>> df = pd.DataFrame({"A": [True, False], "B": [1, 0]})
 |      >>> df
 |             A  B
 |      0   True  1
 |      1  False  0
 |      
 |      >>> df.any(axis='columns')
 |      0    True
 |      1    False
 |      dtype: bool
 |      
 |      Aggregating over the entire DataFrame with ``axis=None``.
 |      
 |      >>> df.any(axis=None)
 |      True
 |      
 |      `any` for an empty DataFrame is an empty Series.
 |      
 |      >>> pd.DataFrame([]).any()
 |      Series([], dtype: bool)
 |  
 |  append(self, other, ignore_index=False, verify_integrity=False, sort=None)
 |      Append rows of `other` to the end of this frame, returning a new
 |      object. Columns not in this frame are added as new columns.
 |      
 |      Parameters
 |      ----------
 |      other : DataFrame or Series/dict-like object, or list of these
 |          The data to append.
 |      ignore_index : boolean, default False
 |          If True, do not use the index labels.
 |      verify_integrity : boolean, default False
 |          If True, raise ValueError on creating index with duplicates.
 |      sort : boolean, default None
 |          Sort columns if the columns of `self` and `other` are not aligned.
 |          The default sorting is deprecated and will change to not-sorting
 |          in a future version of pandas. Explicitly pass ``sort=True`` to
 |          silence the warning and sort. Explicitly pass ``sort=False`` to
 |          silence the warning and not sort.
 |      
 |          .. versionadded:: 0.23.0
 |      
 |      Returns
 |      -------
 |      appended : DataFrame
 |      
 |      Notes
 |      -----
 |      If a list of dict/series is passed and the keys are all contained in
 |      the DataFrame's index, the order of the columns in the resulting
 |      DataFrame will be unchanged.
 |      
 |      Iteratively appending rows to a DataFrame can be more computationally
 |      intensive than a single concatenate. A better solution is to append
 |      those rows to a list and then concatenate the list with the original
 |      DataFrame all at once.
 |      
 |      See also
 |      --------
 |      pandas.concat : General function to concatenate DataFrame, Series
 |          or Panel objects
 |      
 |      Examples
 |      --------
 |      
 |      >>> df = pd.DataFrame([[1, 2], [3, 4]], columns=list('AB'))
 |      >>> df
 |         A  B
 |      0  1  2
 |      1  3  4
 |      >>> df2 = pd.DataFrame([[5, 6], [7, 8]], columns=list('AB'))
 |      >>> df.append(df2)
 |         A  B
 |      0  1  2
 |      1  3  4
 |      0  5  6
 |      1  7  8
 |      
 |      With `ignore_index` set to True:
 |      
 |      >>> df.append(df2, ignore_index=True)
 |         A  B
 |      0  1  2
 |      1  3  4
 |      2  5  6
 |      3  7  8
 |      
 |      The following, while not recommended methods for generating DataFrames,
 |      show two ways to generate a DataFrame from multiple data sources.
 |      
 |      Less efficient:
 |      
 |      >>> df = pd.DataFrame(columns=['A'])
 |      >>> for i in range(5):
 |      ...     df = df.append({'A': i}, ignore_index=True)
 |      >>> df
 |         A
 |      0  0
 |      1  1
 |      2  2
 |      3  3
 |      4  4
 |      
 |      More efficient:
 |      
 |      >>> pd.concat([pd.DataFrame([i], columns=['A']) for i in range(5)],
 |      ...           ignore_index=True)
 |         A
 |      0  0
 |      1  1
 |      2  2
 |      3  3
 |      4  4
 |  
 |  apply(self, func, axis=0, broadcast=None, raw=False, reduce=None, result_type=None, args=(), **kwds)
 |      Apply a function along an axis of the DataFrame.
 |      
 |      Objects passed to the function are Series objects whose index is
 |      either the DataFrame's index (``axis=0``) or the DataFrame's columns
 |      (``axis=1``). By default (``result_type=None``), the final return type
 |      is inferred from the return type of the applied function. Otherwise,
 |      it depends on the `result_type` argument.
 |      
 |      Parameters
 |      ----------
 |      func : function
 |          Function to apply to each column or row.
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          Axis along which the function is applied:
 |      
 |          * 0 or 'index': apply function to each column.
 |          * 1 or 'columns': apply function to each row.
 |      broadcast : bool, optional
 |          Only relevant for aggregation functions:
 |      
 |          * ``False`` or ``None`` : returns a Series whose length is the
 |            length of the index or the number of columns (based on the
 |            `axis` parameter)
 |          * ``True`` : results will be broadcast to the original shape
 |            of the frame, the original index and columns will be retained.
 |      
 |          .. deprecated:: 0.23.0
 |             This argument will be removed in a future version, replaced
 |             by result_type='broadcast'.
 |      
 |      raw : bool, default False
 |          * ``False`` : passes each row or column as a Series to the
 |            function.
 |          * ``True`` : the passed function will receive ndarray objects
 |            instead.
 |            If you are just applying a NumPy reduction function this will
 |            achieve much better performance.
 |      reduce : bool or None, default None
 |          Try to apply reduction procedures. If the DataFrame is empty,
 |          `apply` will use `reduce` to determine whether the result
 |          should be a Series or a DataFrame. If ``reduce=None`` (the
 |          default), `apply`'s return value will be guessed by calling
 |          `func` on an empty Series
 |          (note: while guessing, exceptions raised by `func` will be
 |          ignored).
 |          If ``reduce=True`` a Series will always be returned, and if
 |          ``reduce=False`` a DataFrame will always be returned.
 |      
 |          .. deprecated:: 0.23.0
 |             This argument will be removed in a future version, replaced
 |             by ``result_type='reduce'``.
 |      
 |      result_type : {'expand', 'reduce', 'broadcast', None}, default None
 |          These only act when ``axis=1`` (columns):
 |      
 |          * 'expand' : list-like results will be turned into columns.
 |          * 'reduce' : returns a Series if possible rather than expanding
 |            list-like results. This is the opposite of 'expand'.
 |          * 'broadcast' : results will be broadcast to the original shape
 |            of the DataFrame, the original index and columns will be
 |            retained.
 |      
 |          The default behaviour (None) depends on the return value of the
 |          applied function: list-like results will be returned as a Series
 |          of those. However if the apply function returns a Series these
 |          are expanded to columns.
 |      
 |          .. versionadded:: 0.23.0
 |      
 |      args : tuple
 |          Positional arguments to pass to `func` in addition to the
 |          array/series.
 |      **kwds
 |          Additional keyword arguments to pass as keywords arguments to
 |          `func`.
 |      
 |      Notes
 |      -----
 |      In the current implementation apply calls `func` twice on the
 |      first column/row to decide whether it can take a fast or slow
 |      code path. This can lead to unexpected behavior if `func` has
 |      side-effects, as they will take effect twice for the first
 |      column/row.
 |      
 |      See also
 |      --------
 |      DataFrame.applymap: For elementwise operations
 |      DataFrame.aggregate: only perform aggregating type operations
 |      DataFrame.transform: only perform transformating type operations
 |      
 |      Examples
 |      --------
 |      
 |      >>> df = pd.DataFrame([[4, 9],] * 3, columns=['A', 'B'])
 |      >>> df
 |         A  B
 |      0  4  9
 |      1  4  9
 |      2  4  9
 |      
 |      Using a numpy universal function (in this case the same as
 |      ``np.sqrt(df)``):
 |      
 |      >>> df.apply(np.sqrt)
 |           A    B
 |      0  2.0  3.0
 |      1  2.0  3.0
 |      2  2.0  3.0
 |      
 |      Using a reducing function on either axis
 |      
 |      >>> df.apply(np.sum, axis=0)
 |      A    12
 |      B    27
 |      dtype: int64
 |      
 |      >>> df.apply(np.sum, axis=1)
 |      0    13
 |      1    13
 |      2    13
 |      dtype: int64
 |      
 |      Retuning a list-like will result in a Series
 |      
 |      >>> df.apply(lambda x: [1, 2], axis=1)
 |      0    [1, 2]
 |      1    [1, 2]
 |      2    [1, 2]
 |      dtype: object
 |      
 |      Passing result_type='expand' will expand list-like results
 |      to columns of a Dataframe
 |      
 |      >>> df.apply(lambda x: [1, 2], axis=1, result_type='expand')
 |         0  1
 |      0  1  2
 |      1  1  2
 |      2  1  2
 |      
 |      Returning a Series inside the function is similar to passing
 |      ``result_type='expand'``. The resulting column names
 |      will be the Series index.
 |      
 |      >>> df.apply(lambda x: pd.Series([1, 2], index=['foo', 'bar']), axis=1)
 |         foo  bar
 |      0    1    2
 |      1    1    2
 |      2    1    2
 |      
 |      Passing ``result_type='broadcast'`` will ensure the same shape
 |      result, whether list-like or scalar is returned by the function,
 |      and broadcast it along the axis. The resulting column names will
 |      be the originals.
 |      
 |      >>> df.apply(lambda x: [1, 2], axis=1, result_type='broadcast')
 |         A  B
 |      0  1  2
 |      1  1  2
 |      2  1  2
 |      
 |      Returns
 |      -------
 |      applied : Series or DataFrame
 |  
 |  applymap(self, func)
 |      Apply a function to a Dataframe elementwise.
 |      
 |      This method applies a function that accepts and returns a scalar
 |      to every element of a DataFrame.
 |      
 |      Parameters
 |      ----------
 |      func : callable
 |          Python function, returns a single value from a single value.
 |      
 |      Returns
 |      -------
 |      DataFrame
 |          Transformed DataFrame.
 |      
 |      See also
 |      --------
 |      DataFrame.apply : Apply a function along input axis of DataFrame
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame([[1, 2.12], [3.356, 4.567]])
 |      >>> df
 |             0      1
 |      0  1.000  2.120
 |      1  3.356  4.567
 |      
 |      >>> df.applymap(lambda x: len(str(x)))
 |         0  1
 |      0  3  4
 |      1  5  5
 |      
 |      Note that a vectorized version of `func` often exists, which will
 |      be much faster. You could square each number elementwise.
 |      
 |      >>> df.applymap(lambda x: x**2)
 |                 0          1
 |      0   1.000000   4.494400
 |      1  11.262736  20.857489
 |      
 |      But it's better to avoid applymap in that case.
 |      
 |      >>> df ** 2
 |                 0          1
 |      0   1.000000   4.494400
 |      1  11.262736  20.857489
 |  
 |  assign(self, **kwargs)
 |      Assign new columns to a DataFrame, returning a new object
 |      (a copy) with the new columns added to the original ones.
 |      Existing columns that are re-assigned will be overwritten.
 |      
 |      Parameters
 |      ----------
 |      kwargs : keyword, value pairs
 |          keywords are the column names. If the values are
 |          callable, they are computed on the DataFrame and
 |          assigned to the new columns. The callable must not
 |          change input DataFrame (though pandas doesn't check it).
 |          If the values are not callable, (e.g. a Series, scalar, or array),
 |          they are simply assigned.
 |      
 |      Returns
 |      -------
 |      df : DataFrame
 |          A new DataFrame with the new columns in addition to
 |          all the existing columns.
 |      
 |      Notes
 |      -----
 |      Assigning multiple columns within the same ``assign`` is possible.
 |      For Python 3.6 and above, later items in '\*\*kwargs' may refer to
 |      newly created or modified columns in 'df'; items are computed and
 |      assigned into 'df' in order.  For Python 3.5 and below, the order of
 |      keyword arguments is not specified, you cannot refer to newly created
 |      or modified columns. All items are computed first, and then assigned
 |      in alphabetical order.
 |      
 |      .. versionchanged :: 0.23.0
 |      
 |          Keyword argument order is maintained for Python 3.6 and later.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'A': range(1, 11), 'B': np.random.randn(10)})
 |      
 |      Where the value is a callable, evaluated on `df`:
 |      
 |      >>> df.assign(ln_A = lambda x: np.log(x.A))
 |          A         B      ln_A
 |      0   1  0.426905  0.000000
 |      1   2 -0.780949  0.693147
 |      2   3 -0.418711  1.098612
 |      3   4 -0.269708  1.386294
 |      4   5 -0.274002  1.609438
 |      5   6 -0.500792  1.791759
 |      6   7  1.649697  1.945910
 |      7   8 -1.495604  2.079442
 |      8   9  0.549296  2.197225
 |      9  10 -0.758542  2.302585
 |      
 |      Where the value already exists and is inserted:
 |      
 |      >>> newcol = np.log(df['A'])
 |      >>> df.assign(ln_A=newcol)
 |          A         B      ln_A
 |      0   1  0.426905  0.000000
 |      1   2 -0.780949  0.693147
 |      2   3 -0.418711  1.098612
 |      3   4 -0.269708  1.386294
 |      4   5 -0.274002  1.609438
 |      5   6 -0.500792  1.791759
 |      6   7  1.649697  1.945910
 |      7   8 -1.495604  2.079442
 |      8   9  0.549296  2.197225
 |      9  10 -0.758542  2.302585
 |      
 |      Where the keyword arguments depend on each other
 |      
 |      >>> df = pd.DataFrame({'A': [1, 2, 3]})
 |      
 |      >>> df.assign(B=df.A, C=lambda x:x['A']+ x['B'])
 |          A  B  C
 |       0  1  1  2
 |       1  2  2  4
 |       2  3  3  6
 |  
 |  boxplot = boxplot_frame(self, column=None, by=None, ax=None, fontsize=None, rot=0, grid=True, figsize=None, layout=None, return_type=None, **kwds)
 |      Make a box plot from DataFrame columns.
 |      
 |      Make a box-and-whisker plot from DataFrame columns, optionally grouped
 |      by some other columns. A box plot is a method for graphically depicting
 |      groups of numerical data through their quartiles.
 |      The box extends from the Q1 to Q3 quartile values of the data,
 |      with a line at the median (Q2). The whiskers extend from the edges
 |      of box to show the range of the data. The position of the whiskers
 |      is set by default to `1.5 * IQR (IQR = Q3 - Q1)` from the edges of the box.
 |      Outlier points are those past the end of the whiskers.
 |      
 |      For further details see
 |      Wikipedia's entry for `boxplot <https://en.wikipedia.org/wiki/Box_plot>`_.
 |      
 |      Parameters
 |      ----------
 |      column : str or list of str, optional
 |          Column name or list of names, or vector.
 |          Can be any valid input to :meth:`pandas.DataFrame.groupby`.
 |      by : str or array-like, optional
 |          Column in the DataFrame to :meth:`pandas.DataFrame.groupby`.
 |          One box-plot will be done per value of columns in `by`.
 |      ax : object of class matplotlib.axes.Axes, optional
 |          The matplotlib axes to be used by boxplot.
 |      fontsize : float or str
 |          Tick label font size in points or as a string (e.g., `large`).
 |      rot : int or float, default 0
 |          The rotation angle of labels (in degrees)
 |          with respect to the screen coordinate sytem.
 |      grid : boolean, default True
 |          Setting this to True will show the grid.
 |      figsize : A tuple (width, height) in inches
 |          The size of the figure to create in matplotlib.
 |      layout : tuple (rows, columns), optional
 |          For example, (3, 5) will display the subplots
 |          using 3 columns and 5 rows, starting from the top-left.
 |      return_type : {'axes', 'dict', 'both'} or None, default 'axes'
 |          The kind of object to return. The default is ``axes``.
 |      
 |          * 'axes' returns the matplotlib axes the boxplot is drawn on.
 |          * 'dict' returns a dictionary whose values are the matplotlib
 |            Lines of the boxplot.
 |          * 'both' returns a namedtuple with the axes and dict.
 |          * when grouping with ``by``, a Series mapping columns to
 |            ``return_type`` is returned.
 |      
 |            If ``return_type`` is `None`, a NumPy array
 |            of axes with the same shape as ``layout`` is returned.
 |      **kwds
 |          All other plotting keyword arguments to be passed to
 |          :func:`matplotlib.pyplot.boxplot`.
 |      
 |      Returns
 |      -------
 |      result :
 |      
 |          The return type depends on the `return_type` parameter:
 |      
 |          * 'axes' : object of class matplotlib.axes.Axes
 |          * 'dict' : dict of matplotlib.lines.Line2D objects
 |          * 'both' : a nametuple with strucure (ax, lines)
 |      
 |          For data grouped with ``by``:
 |      
 |          * :class:`~pandas.Series`
 |          * :class:`~numpy.array` (for ``return_type = None``)
 |      
 |      See Also
 |      --------
 |      Series.plot.hist: Make a histogram.
 |      matplotlib.pyplot.boxplot : Matplotlib equivalent plot.
 |      
 |      Notes
 |      -----
 |      Use ``return_type='dict'`` when you want to tweak the appearance
 |      of the lines after plotting. In this case a dict containing the Lines
 |      making up the boxes, caps, fliers, medians, and whiskers is returned.
 |      
 |      Examples
 |      --------
 |      
 |      Boxplots can be created for every column in the dataframe
 |      by ``df.boxplot()`` or indicating the columns to be used:
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> np.random.seed(1234)
 |          >>> df = pd.DataFrame(np.random.randn(10,4),
 |          ...                   columns=['Col1', 'Col2', 'Col3', 'Col4'])
 |          >>> boxplot = df.boxplot(column=['Col1', 'Col2', 'Col3'])
 |      
 |      Boxplots of variables distributions grouped by the values of a third
 |      variable can be created using the option ``by``. For instance:
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> df = pd.DataFrame(np.random.randn(10, 2),
 |          ...                   columns=['Col1', 'Col2'])
 |          >>> df['X'] = pd.Series(['A', 'A', 'A', 'A', 'A',
 |          ...                      'B', 'B', 'B', 'B', 'B'])
 |          >>> boxplot = df.boxplot(by='X')
 |      
 |      A list of strings (i.e. ``['X', 'Y']``) can be passed to boxplot
 |      in order to group the data by combination of the variables in the x-axis:
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> df = pd.DataFrame(np.random.randn(10,3),
 |          ...                   columns=['Col1', 'Col2', 'Col3'])
 |          >>> df['X'] = pd.Series(['A', 'A', 'A', 'A', 'A',
 |          ...                      'B', 'B', 'B', 'B', 'B'])
 |          >>> df['Y'] = pd.Series(['A', 'B', 'A', 'B', 'A',
 |          ...                      'B', 'A', 'B', 'A', 'B'])
 |          >>> boxplot = df.boxplot(column=['Col1', 'Col2'], by=['X', 'Y'])
 |      
 |      The layout of boxplot can be adjusted giving a tuple to ``layout``:
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> boxplot = df.boxplot(column=['Col1', 'Col2'], by='X',
 |          ...                      layout=(2, 1))
 |      
 |      Additional formatting can be done to the boxplot, like suppressing the grid
 |      (``grid=False``), rotating the labels in the x-axis (i.e. ``rot=45``)
 |      or changing the fontsize (i.e. ``fontsize=15``):
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          >>> boxplot = df.boxplot(grid=False, rot=45, fontsize=15)
 |      
 |      The parameter ``return_type`` can be used to select the type of element
 |      returned by `boxplot`.  When ``return_type='axes'`` is selected,
 |      the matplotlib axes on which the boxplot is drawn are returned:
 |      
 |          >>> boxplot = df.boxplot(column=['Col1','Col2'], return_type='axes')
 |          >>> type(boxplot)
 |          <class 'matplotlib.axes._subplots.AxesSubplot'>
 |      
 |      When grouping with ``by``, a Series mapping columns to ``return_type``
 |      is returned:
 |      
 |          >>> boxplot = df.boxplot(column=['Col1', 'Col2'], by='X',
 |          ...                      return_type='axes')
 |          >>> type(boxplot)
 |          <class 'pandas.core.series.Series'>
 |      
 |      If ``return_type`` is `None`, a NumPy array of axes with the same shape
 |      as ``layout`` is returned:
 |      
 |          >>> boxplot =  df.boxplot(column=['Col1', 'Col2'], by='X',
 |          ...                       return_type=None)
 |          >>> type(boxplot)
 |          <class 'numpy.ndarray'>
 |  
 |  combine(self, other, func, fill_value=None, overwrite=True)
 |      Add two DataFrame objects and do not propagate NaN values, so if for a
 |      (column, time) one frame is missing a value, it will default to the
 |      other frame's value (which might be NaN as well)
 |      
 |      Parameters
 |      ----------
 |      other : DataFrame
 |      func : function
 |          Function that takes two series as inputs and return a Series or a
 |          scalar
 |      fill_value : scalar value
 |      overwrite : boolean, default True
 |          If True then overwrite values for common keys in the calling frame
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      >>> df1 = DataFrame({'A': [0, 0], 'B': [4, 4]})
 |      >>> df2 = DataFrame({'A': [1, 1], 'B': [3, 3]})
 |      >>> df1.combine(df2, lambda s1, s2: s1 if s1.sum() < s2.sum() else s2)
 |         A  B
 |      0  0  3
 |      1  0  3
 |      
 |      See Also
 |      --------
 |      DataFrame.combine_first : Combine two DataFrame objects and default to
 |          non-null values in frame calling the method
 |  
 |  combine_first(self, other)
 |      Combine two DataFrame objects and default to non-null values in frame
 |      calling the method. Result index columns will be the union of the
 |      respective indexes and columns
 |      
 |      Parameters
 |      ----------
 |      other : DataFrame
 |      
 |      Returns
 |      -------
 |      combined : DataFrame
 |      
 |      Examples
 |      --------
 |      df1's values prioritized, use values from df2 to fill holes:
 |      
 |      >>> df1 = pd.DataFrame([[1, np.nan]])
 |      >>> df2 = pd.DataFrame([[3, 4]])
 |      >>> df1.combine_first(df2)
 |         0    1
 |      0  1  4.0
 |      
 |      See Also
 |      --------
 |      DataFrame.combine : Perform series-wise operation on two DataFrames
 |          using a given function
 |  
 |  compound(self, axis=None, skipna=None, level=None)
 |      Return the compound percentage of the values for the requested axis
 |      
 |      Parameters
 |      ----------
 |      axis : {index (0), columns (1)}
 |      skipna : boolean, default True
 |          Exclude NA/null values when computing the result.
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series
 |      numeric_only : boolean, default None
 |          Include only float, int, boolean columns. If None, will attempt to use
 |          everything, then use only numeric data. Not implemented for Series.
 |      
 |      Returns
 |      -------
 |      compounded : Series or DataFrame (if level specified)
 |  
 |  corr(self, method='pearson', min_periods=1)
 |      Compute pairwise correlation of columns, excluding NA/null values
 |      
 |      Parameters
 |      ----------
 |      method : {'pearson', 'kendall', 'spearman'}
 |          * pearson : standard correlation coefficient
 |          * kendall : Kendall Tau correlation coefficient
 |          * spearman : Spearman rank correlation
 |      min_periods : int, optional
 |          Minimum number of observations required per pair of columns
 |          to have a valid result. Currently only available for pearson
 |          and spearman correlation
 |      
 |      Returns
 |      -------
 |      y : DataFrame
 |  
 |  corrwith(self, other, axis=0, drop=False)
 |      Compute pairwise correlation between rows or columns of two DataFrame
 |      objects.
 |      
 |      Parameters
 |      ----------
 |      other : DataFrame, Series
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          0 or 'index' to compute column-wise, 1 or 'columns' for row-wise
 |      drop : boolean, default False
 |          Drop missing indices from result, default returns union of all
 |      
 |      Returns
 |      -------
 |      correls : Series
 |  
 |  count(self, axis=0, level=None, numeric_only=False)
 |      Count non-NA cells for each column or row.
 |      
 |      The values `None`, `NaN`, `NaT`, and optionally `numpy.inf` (depending
 |      on `pandas.options.mode.use_inf_as_na`) are considered NA.
 |      
 |      Parameters
 |      ----------
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          If 0 or 'index' counts are generated for each column.
 |          If 1 or 'columns' counts are generated for each **row**.
 |      level : int or str, optional
 |          If the axis is a `MultiIndex` (hierarchical), count along a
 |          particular `level`, collapsing into a `DataFrame`.
 |          A `str` specifies the level name.
 |      numeric_only : boolean, default False
 |          Include only `float`, `int` or `boolean` data.
 |      
 |      Returns
 |      -------
 |      Series or DataFrame
 |          For each column/row the number of non-NA/null entries.
 |          If `level` is specified returns a `DataFrame`.
 |      
 |      See Also
 |      --------
 |      Series.count: number of non-NA elements in a Series
 |      DataFrame.shape: number of DataFrame rows and columns (including NA
 |          elements)
 |      DataFrame.isna: boolean same-sized DataFrame showing places of NA
 |          elements
 |      
 |      Examples
 |      --------
 |      Constructing DataFrame from a dictionary:
 |      
 |      >>> df = pd.DataFrame({"Person":
 |      ...                    ["John", "Myla", None, "John", "Myla"],
 |      ...                    "Age": [24., np.nan, 21., 33, 26],
 |      ...                    "Single": [False, True, True, True, False]})
 |      >>> df
 |         Person   Age  Single
 |      0    John  24.0   False
 |      1    Myla   NaN    True
 |      2    None  21.0    True
 |      3    John  33.0    True
 |      4    Myla  26.0   False
 |      
 |      Notice the uncounted NA values:
 |      
 |      >>> df.count()
 |      Person    4
 |      Age       4
 |      Single    5
 |      dtype: int64
 |      
 |      Counts for each **row**:
 |      
 |      >>> df.count(axis='columns')
 |      0    3
 |      1    2
 |      2    2
 |      3    3
 |      4    3
 |      dtype: int64
 |      
 |      Counts for one level of a `MultiIndex`:
 |      
 |      >>> df.set_index(["Person", "Single"]).count(level="Person")
 |              Age
 |      Person
 |      John      2
 |      Myla      1
 |  
 |  cov(self, min_periods=None)
 |      Compute pairwise covariance of columns, excluding NA/null values.
 |      
 |      Compute the pairwise covariance among the series of a DataFrame.
 |      The returned data frame is the `covariance matrix
 |      <https://en.wikipedia.org/wiki/Covariance_matrix>`__ of the columns
 |      of the DataFrame.
 |      
 |      Both NA and null values are automatically excluded from the
 |      calculation. (See the note below about bias from missing values.)
 |      A threshold can be set for the minimum number of
 |      observations for each value created. Comparisons with observations
 |      below this threshold will be returned as ``NaN``.
 |      
 |      This method is generally used for the analysis of time series data to
 |      understand the relationship between different measures
 |      across time.
 |      
 |      Parameters
 |      ----------
 |      min_periods : int, optional
 |          Minimum number of observations required per pair of columns
 |          to have a valid result.
 |      
 |      Returns
 |      -------
 |      DataFrame
 |          The covariance matrix of the series of the DataFrame.
 |      
 |      See Also
 |      --------
 |      pandas.Series.cov : compute covariance with another Series
 |      pandas.core.window.EWM.cov: expoential weighted sample covariance
 |      pandas.core.window.Expanding.cov : expanding sample covariance
 |      pandas.core.window.Rolling.cov : rolling sample covariance
 |      
 |      Notes
 |      -----
 |      Returns the covariance matrix of the DataFrame's time series.
 |      The covariance is normalized by N-1.
 |      
 |      For DataFrames that have Series that are missing data (assuming that
 |      data is `missing at random
 |      <https://en.wikipedia.org/wiki/Missing_data#Missing_at_random>`__)
 |      the returned covariance matrix will be an unbiased estimate
 |      of the variance and covariance between the member Series.
 |      
 |      However, for many applications this estimate may not be acceptable
 |      because the estimate covariance matrix is not guaranteed to be positive
 |      semi-definite. This could lead to estimate correlations having
 |      absolute values which are greater than one, and/or a non-invertible
 |      covariance matrix. See `Estimation of covariance matrices
 |      <http://en.wikipedia.org/w/index.php?title=Estimation_of_covariance_
 |      matrices>`__ for more details.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame([(1, 2), (0, 3), (2, 0), (1, 1)],
 |      ...                   columns=['dogs', 'cats'])
 |      >>> df.cov()
 |                dogs      cats
 |      dogs  0.666667 -1.000000
 |      cats -1.000000  1.666667
 |      
 |      >>> np.random.seed(42)
 |      >>> df = pd.DataFrame(np.random.randn(1000, 5),
 |      ...                   columns=['a', 'b', 'c', 'd', 'e'])
 |      >>> df.cov()
 |                a         b         c         d         e
 |      a  0.998438 -0.020161  0.059277 -0.008943  0.014144
 |      b -0.020161  1.059352 -0.008543 -0.024738  0.009826
 |      c  0.059277 -0.008543  1.010670 -0.001486 -0.000271
 |      d -0.008943 -0.024738 -0.001486  0.921297 -0.013692
 |      e  0.014144  0.009826 -0.000271 -0.013692  0.977795
 |      
 |      **Minimum number of periods**
 |      
 |      This method also supports an optional ``min_periods`` keyword
 |      that specifies the required minimum number of non-NA observations for
 |      each column pair in order to have a valid result:
 |      
 |      >>> np.random.seed(42)
 |      >>> df = pd.DataFrame(np.random.randn(20, 3),
 |      ...                   columns=['a', 'b', 'c'])
 |      >>> df.loc[df.index[:5], 'a'] = np.nan
 |      >>> df.loc[df.index[5:10], 'b'] = np.nan
 |      >>> df.cov(min_periods=12)
 |                a         b         c
 |      a  0.316741       NaN -0.150812
 |      b       NaN  1.248003  0.191417
 |      c -0.150812  0.191417  0.895202
 |  
 |  cummax(self, axis=None, skipna=True, *args, **kwargs)
 |      Return cumulative maximum over a DataFrame or Series axis.
 |      
 |      Returns a DataFrame or Series of the same size containing the cumulative
 |      maximum.
 |      
 |      Parameters
 |      ----------
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          The index or the name of the axis. 0 is equivalent to None or 'index'.
 |      skipna : boolean, default True
 |          Exclude NA/null values. If an entire row/column is NA, the result
 |          will be NA.
 |      *args, **kwargs :
 |          Additional keywords have no effect but might be accepted for
 |          compatibility with NumPy.
 |      
 |      Returns
 |      -------
 |      cummax : Series or DataFrame
 |      
 |      Examples
 |      --------
 |      **Series**
 |      
 |      >>> s = pd.Series([2, np.nan, 5, -1, 0])
 |      >>> s
 |      0    2.0
 |      1    NaN
 |      2    5.0
 |      3   -1.0
 |      4    0.0
 |      dtype: float64
 |      
 |      By default, NA values are ignored.
 |      
 |      >>> s.cummax()
 |      0    2.0
 |      1    NaN
 |      2    5.0
 |      3    5.0
 |      4    5.0
 |      dtype: float64
 |      
 |      To include NA values in the operation, use ``skipna=False``
 |      
 |      >>> s.cummax(skipna=False)
 |      0    2.0
 |      1    NaN
 |      2    NaN
 |      3    NaN
 |      4    NaN
 |      dtype: float64
 |      
 |      **DataFrame**
 |      
 |      >>> df = pd.DataFrame([[2.0, 1.0],
 |      ...                    [3.0, np.nan],
 |      ...                    [1.0, 0.0]],
 |      ...                    columns=list('AB'))
 |      >>> df
 |           A    B
 |      0  2.0  1.0
 |      1  3.0  NaN
 |      2  1.0  0.0
 |      
 |      By default, iterates over rows and finds the maximum
 |      in each column. This is equivalent to ``axis=None`` or ``axis='index'``.
 |      
 |      >>> df.cummax()
 |           A    B
 |      0  2.0  1.0
 |      1  3.0  NaN
 |      2  3.0  1.0
 |      
 |      To iterate over columns and find the maximum in each row,
 |      use ``axis=1``
 |      
 |      >>> df.cummax(axis=1)
 |           A    B
 |      0  2.0  2.0
 |      1  3.0  NaN
 |      2  1.0  1.0
 |      
 |      See also
 |      --------
 |      pandas.core.window.Expanding.max : Similar functionality
 |          but ignores ``NaN`` values.
 |      DataFrame.max : Return the maximum over
 |          DataFrame axis.
 |      DataFrame.cummax : Return cumulative maximum over DataFrame axis.
 |      DataFrame.cummin : Return cumulative minimum over DataFrame axis.
 |      DataFrame.cumsum : Return cumulative sum over DataFrame axis.
 |      DataFrame.cumprod : Return cumulative product over DataFrame axis.
 |  
 |  cummin(self, axis=None, skipna=True, *args, **kwargs)
 |      Return cumulative minimum over a DataFrame or Series axis.
 |      
 |      Returns a DataFrame or Series of the same size containing the cumulative
 |      minimum.
 |      
 |      Parameters
 |      ----------
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          The index or the name of the axis. 0 is equivalent to None or 'index'.
 |      skipna : boolean, default True
 |          Exclude NA/null values. If an entire row/column is NA, the result
 |          will be NA.
 |      *args, **kwargs :
 |          Additional keywords have no effect but might be accepted for
 |          compatibility with NumPy.
 |      
 |      Returns
 |      -------
 |      cummin : Series or DataFrame
 |      
 |      Examples
 |      --------
 |      **Series**
 |      
 |      >>> s = pd.Series([2, np.nan, 5, -1, 0])
 |      >>> s
 |      0    2.0
 |      1    NaN
 |      2    5.0
 |      3   -1.0
 |      4    0.0
 |      dtype: float64
 |      
 |      By default, NA values are ignored.
 |      
 |      >>> s.cummin()
 |      0    2.0
 |      1    NaN
 |      2    2.0
 |      3   -1.0
 |      4   -1.0
 |      dtype: float64
 |      
 |      To include NA values in the operation, use ``skipna=False``
 |      
 |      >>> s.cummin(skipna=False)
 |      0    2.0
 |      1    NaN
 |      2    NaN
 |      3    NaN
 |      4    NaN
 |      dtype: float64
 |      
 |      **DataFrame**
 |      
 |      >>> df = pd.DataFrame([[2.0, 1.0],
 |      ...                    [3.0, np.nan],
 |      ...                    [1.0, 0.0]],
 |      ...                    columns=list('AB'))
 |      >>> df
 |           A    B
 |      0  2.0  1.0
 |      1  3.0  NaN
 |      2  1.0  0.0
 |      
 |      By default, iterates over rows and finds the minimum
 |      in each column. This is equivalent to ``axis=None`` or ``axis='index'``.
 |      
 |      >>> df.cummin()
 |           A    B
 |      0  2.0  1.0
 |      1  2.0  NaN
 |      2  1.0  0.0
 |      
 |      To iterate over columns and find the minimum in each row,
 |      use ``axis=1``
 |      
 |      >>> df.cummin(axis=1)
 |           A    B
 |      0  2.0  1.0
 |      1  3.0  NaN
 |      2  1.0  0.0
 |      
 |      See also
 |      --------
 |      pandas.core.window.Expanding.min : Similar functionality
 |          but ignores ``NaN`` values.
 |      DataFrame.min : Return the minimum over
 |          DataFrame axis.
 |      DataFrame.cummax : Return cumulative maximum over DataFrame axis.
 |      DataFrame.cummin : Return cumulative minimum over DataFrame axis.
 |      DataFrame.cumsum : Return cumulative sum over DataFrame axis.
 |      DataFrame.cumprod : Return cumulative product over DataFrame axis.
 |  
 |  cumprod(self, axis=None, skipna=True, *args, **kwargs)
 |      Return cumulative product over a DataFrame or Series axis.
 |      
 |      Returns a DataFrame or Series of the same size containing the cumulative
 |      product.
 |      
 |      Parameters
 |      ----------
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          The index or the name of the axis. 0 is equivalent to None or 'index'.
 |      skipna : boolean, default True
 |          Exclude NA/null values. If an entire row/column is NA, the result
 |          will be NA.
 |      *args, **kwargs :
 |          Additional keywords have no effect but might be accepted for
 |          compatibility with NumPy.
 |      
 |      Returns
 |      -------
 |      cumprod : Series or DataFrame
 |      
 |      Examples
 |      --------
 |      **Series**
 |      
 |      >>> s = pd.Series([2, np.nan, 5, -1, 0])
 |      >>> s
 |      0    2.0
 |      1    NaN
 |      2    5.0
 |      3   -1.0
 |      4    0.0
 |      dtype: float64
 |      
 |      By default, NA values are ignored.
 |      
 |      >>> s.cumprod()
 |      0     2.0
 |      1     NaN
 |      2    10.0
 |      3   -10.0
 |      4    -0.0
 |      dtype: float64
 |      
 |      To include NA values in the operation, use ``skipna=False``
 |      
 |      >>> s.cumprod(skipna=False)
 |      0    2.0
 |      1    NaN
 |      2    NaN
 |      3    NaN
 |      4    NaN
 |      dtype: float64
 |      
 |      **DataFrame**
 |      
 |      >>> df = pd.DataFrame([[2.0, 1.0],
 |      ...                    [3.0, np.nan],
 |      ...                    [1.0, 0.0]],
 |      ...                    columns=list('AB'))
 |      >>> df
 |           A    B
 |      0  2.0  1.0
 |      1  3.0  NaN
 |      2  1.0  0.0
 |      
 |      By default, iterates over rows and finds the product
 |      in each column. This is equivalent to ``axis=None`` or ``axis='index'``.
 |      
 |      >>> df.cumprod()
 |           A    B
 |      0  2.0  1.0
 |      1  6.0  NaN
 |      2  6.0  0.0
 |      
 |      To iterate over columns and find the product in each row,
 |      use ``axis=1``
 |      
 |      >>> df.cumprod(axis=1)
 |           A    B
 |      0  2.0  2.0
 |      1  3.0  NaN
 |      2  1.0  0.0
 |      
 |      See also
 |      --------
 |      pandas.core.window.Expanding.prod : Similar functionality
 |          but ignores ``NaN`` values.
 |      DataFrame.prod : Return the product over
 |          DataFrame axis.
 |      DataFrame.cummax : Return cumulative maximum over DataFrame axis.
 |      DataFrame.cummin : Return cumulative minimum over DataFrame axis.
 |      DataFrame.cumsum : Return cumulative sum over DataFrame axis.
 |      DataFrame.cumprod : Return cumulative product over DataFrame axis.
 |  
 |  cumsum(self, axis=None, skipna=True, *args, **kwargs)
 |      Return cumulative sum over a DataFrame or Series axis.
 |      
 |      Returns a DataFrame or Series of the same size containing the cumulative
 |      sum.
 |      
 |      Parameters
 |      ----------
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          The index or the name of the axis. 0 is equivalent to None or 'index'.
 |      skipna : boolean, default True
 |          Exclude NA/null values. If an entire row/column is NA, the result
 |          will be NA.
 |      *args, **kwargs :
 |          Additional keywords have no effect but might be accepted for
 |          compatibility with NumPy.
 |      
 |      Returns
 |      -------
 |      cumsum : Series or DataFrame
 |      
 |      Examples
 |      --------
 |      **Series**
 |      
 |      >>> s = pd.Series([2, np.nan, 5, -1, 0])
 |      >>> s
 |      0    2.0
 |      1    NaN
 |      2    5.0
 |      3   -1.0
 |      4    0.0
 |      dtype: float64
 |      
 |      By default, NA values are ignored.
 |      
 |      >>> s.cumsum()
 |      0    2.0
 |      1    NaN
 |      2    7.0
 |      3    6.0
 |      4    6.0
 |      dtype: float64
 |      
 |      To include NA values in the operation, use ``skipna=False``
 |      
 |      >>> s.cumsum(skipna=False)
 |      0    2.0
 |      1    NaN
 |      2    NaN
 |      3    NaN
 |      4    NaN
 |      dtype: float64
 |      
 |      **DataFrame**
 |      
 |      >>> df = pd.DataFrame([[2.0, 1.0],
 |      ...                    [3.0, np.nan],
 |      ...                    [1.0, 0.0]],
 |      ...                    columns=list('AB'))
 |      >>> df
 |           A    B
 |      0  2.0  1.0
 |      1  3.0  NaN
 |      2  1.0  0.0
 |      
 |      By default, iterates over rows and finds the sum
 |      in each column. This is equivalent to ``axis=None`` or ``axis='index'``.
 |      
 |      >>> df.cumsum()
 |           A    B
 |      0  2.0  1.0
 |      1  5.0  NaN
 |      2  6.0  1.0
 |      
 |      To iterate over columns and find the sum in each row,
 |      use ``axis=1``
 |      
 |      >>> df.cumsum(axis=1)
 |           A    B
 |      0  2.0  3.0
 |      1  3.0  NaN
 |      2  1.0  1.0
 |      
 |      See also
 |      --------
 |      pandas.core.window.Expanding.sum : Similar functionality
 |          but ignores ``NaN`` values.
 |      DataFrame.sum : Return the sum over
 |          DataFrame axis.
 |      DataFrame.cummax : Return cumulative maximum over DataFrame axis.
 |      DataFrame.cummin : Return cumulative minimum over DataFrame axis.
 |      DataFrame.cumsum : Return cumulative sum over DataFrame axis.
 |      DataFrame.cumprod : Return cumulative product over DataFrame axis.
 |  
 |  diff(self, periods=1, axis=0)
 |      First discrete difference of element.
 |      
 |      Calculates the difference of a DataFrame element compared with another
 |      element in the DataFrame (default is the element in the same column
 |      of the previous row).
 |      
 |      Parameters
 |      ----------
 |      periods : int, default 1
 |          Periods to shift for calculating difference, accepts negative
 |          values.
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          Take difference over rows (0) or columns (1).
 |      
 |          .. versionadded:: 0.16.1.
 |      
 |      Returns
 |      -------
 |      diffed : DataFrame
 |      
 |      See Also
 |      --------
 |      Series.diff: First discrete difference for a Series.
 |      DataFrame.pct_change: Percent change over given number of periods.
 |      DataFrame.shift: Shift index by desired number of periods with an
 |          optional time freq.
 |      
 |      Examples
 |      --------
 |      Difference with previous row
 |      
 |      >>> df = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6],
 |      ...                    'b': [1, 1, 2, 3, 5, 8],
 |      ...                    'c': [1, 4, 9, 16, 25, 36]})
 |      >>> df
 |         a  b   c
 |      0  1  1   1
 |      1  2  1   4
 |      2  3  2   9
 |      3  4  3  16
 |      4  5  5  25
 |      5  6  8  36
 |      
 |      >>> df.diff()
 |           a    b     c
 |      0  NaN  NaN   NaN
 |      1  1.0  0.0   3.0
 |      2  1.0  1.0   5.0
 |      3  1.0  1.0   7.0
 |      4  1.0  2.0   9.0
 |      5  1.0  3.0  11.0
 |      
 |      Difference with previous column
 |      
 |      >>> df.diff(axis=1)
 |          a    b     c
 |      0 NaN  0.0   0.0
 |      1 NaN -1.0   3.0
 |      2 NaN -1.0   7.0
 |      3 NaN -1.0  13.0
 |      4 NaN  0.0  20.0
 |      5 NaN  2.0  28.0
 |      
 |      Difference with 3rd previous row
 |      
 |      >>> df.diff(periods=3)
 |           a    b     c
 |      0  NaN  NaN   NaN
 |      1  NaN  NaN   NaN
 |      2  NaN  NaN   NaN
 |      3  3.0  2.0  15.0
 |      4  3.0  4.0  21.0
 |      5  3.0  6.0  27.0
 |      
 |      Difference with following row
 |      
 |      >>> df.diff(periods=-1)
 |           a    b     c
 |      0 -1.0  0.0  -3.0
 |      1 -1.0 -1.0  -5.0
 |      2 -1.0 -1.0  -7.0
 |      3 -1.0 -2.0  -9.0
 |      4 -1.0 -3.0 -11.0
 |      5  NaN  NaN   NaN
 |  
 |  div = truediv(self, other, axis='columns', level=None, fill_value=None)
 |  
 |  divide = truediv(self, other, axis='columns', level=None, fill_value=None)
 |  
 |  dot(self, other)
 |      Matrix multiplication with DataFrame or Series objects.  Can also be
 |      called using `self @ other` in Python >= 3.5.
 |      
 |      Parameters
 |      ----------
 |      other : DataFrame or Series
 |      
 |      Returns
 |      -------
 |      dot_product : DataFrame or Series
 |  
 |  drop(self, labels=None, axis=0, index=None, columns=None, level=None, inplace=False, errors='raise')
 |      Drop specified labels from rows or columns.
 |      
 |      Remove rows or columns by specifying label names and corresponding
 |      axis, or by specifying directly index or column names. When using a
 |      multi-index, labels on different levels can be removed by specifying
 |      the level.
 |      
 |      Parameters
 |      ----------
 |      labels : single label or list-like
 |          Index or column labels to drop.
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          Whether to drop labels from the index (0 or 'index') or
 |          columns (1 or 'columns').
 |      index, columns : single label or list-like
 |          Alternative to specifying axis (``labels, axis=1``
 |          is equivalent to ``columns=labels``).
 |      
 |          .. versionadded:: 0.21.0
 |      level : int or level name, optional
 |          For MultiIndex, level from which the labels will be removed.
 |      inplace : bool, default False
 |          If True, do operation inplace and return None.
 |      errors : {'ignore', 'raise'}, default 'raise'
 |          If 'ignore', suppress error and only existing labels are
 |          dropped.
 |      
 |      Returns
 |      -------
 |      dropped : pandas.DataFrame
 |      
 |      See Also
 |      --------
 |      DataFrame.loc : Label-location based indexer for selection by label.
 |      DataFrame.dropna : Return DataFrame with labels on given axis omitted
 |          where (all or any) data are missing
 |      DataFrame.drop_duplicates : Return DataFrame with duplicate rows
 |          removed, optionally only considering certain columns
 |      Series.drop : Return Series with specified index labels removed.
 |      
 |      Raises
 |      ------
 |      KeyError
 |          If none of the labels are found in the selected axis
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame(np.arange(12).reshape(3,4),
 |      ...                   columns=['A', 'B', 'C', 'D'])
 |      >>> df
 |         A  B   C   D
 |      0  0  1   2   3
 |      1  4  5   6   7
 |      2  8  9  10  11
 |      
 |      Drop columns
 |      
 |      >>> df.drop(['B', 'C'], axis=1)
 |         A   D
 |      0  0   3
 |      1  4   7
 |      2  8  11
 |      
 |      >>> df.drop(columns=['B', 'C'])
 |         A   D
 |      0  0   3
 |      1  4   7
 |      2  8  11
 |      
 |      Drop a row by index
 |      
 |      >>> df.drop([0, 1])
 |         A  B   C   D
 |      2  8  9  10  11
 |      
 |      Drop columns and/or rows of MultiIndex DataFrame
 |      
 |      >>> midx = pd.MultiIndex(levels=[['lama', 'cow', 'falcon'],
 |      ...                              ['speed', 'weight', 'length']],
 |      ...                      labels=[[0, 0, 0, 1, 1, 1, 2, 2, 2],
 |      ...                              [0, 1, 2, 0, 1, 2, 0, 1, 2]])
 |      >>> df = pd.DataFrame(index=midx, columns=['big', 'small'],
 |      ...                   data=[[45, 30], [200, 100], [1.5, 1], [30, 20],
 |      ...                         [250, 150], [1.5, 0.8], [320, 250],
 |      ...                         [1, 0.8], [0.3,0.2]])
 |      >>> df
 |                      big     small
 |      lama    speed   45.0    30.0
 |              weight  200.0   100.0
 |              length  1.5     1.0
 |      cow     speed   30.0    20.0
 |              weight  250.0   150.0
 |              length  1.5     0.8
 |      falcon  speed   320.0   250.0
 |              weight  1.0     0.8
 |              length  0.3     0.2
 |      
 |      >>> df.drop(index='cow', columns='small')
 |                      big
 |      lama    speed   45.0
 |              weight  200.0
 |              length  1.5
 |      falcon  speed   320.0
 |              weight  1.0
 |              length  0.3
 |      
 |      >>> df.drop(index='length', level=1)
 |                      big     small
 |      lama    speed   45.0    30.0
 |              weight  200.0   100.0
 |      cow     speed   30.0    20.0
 |              weight  250.0   150.0
 |      falcon  speed   320.0   250.0
 |              weight  1.0     0.8
 |  
 |  drop_duplicates(self, subset=None, keep='first', inplace=False)
 |      Return DataFrame with duplicate rows removed, optionally only
 |      considering certain columns
 |      
 |      Parameters
 |      ----------
 |      subset : column label or sequence of labels, optional
 |          Only consider certain columns for identifying duplicates, by
 |          default use all of the columns
 |      keep : {'first', 'last', False}, default 'first'
 |          - ``first`` : Drop duplicates except for the first occurrence.
 |          - ``last`` : Drop duplicates except for the last occurrence.
 |          - False : Drop all duplicates.
 |      inplace : boolean, default False
 |          Whether to drop duplicates in place or to return a copy
 |      
 |      Returns
 |      -------
 |      deduplicated : DataFrame
 |  
 |  dropna(self, axis=0, how='any', thresh=None, subset=None, inplace=False)
 |      Remove missing values.
 |      
 |      See the :ref:`User Guide <missing_data>` for more on which values are
 |      considered missing, and how to work with missing data.
 |      
 |      Parameters
 |      ----------
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          Determine if rows or columns which contain missing values are
 |          removed.
 |      
 |          * 0, or 'index' : Drop rows which contain missing values.
 |          * 1, or 'columns' : Drop columns which contain missing value.
 |      
 |          .. deprecated:: 0.23.0: Pass tuple or list to drop on multiple
 |          axes.
 |      how : {'any', 'all'}, default 'any'
 |          Determine if row or column is removed from DataFrame, when we have
 |          at least one NA or all NA.
 |      
 |          * 'any' : If any NA values are present, drop that row or column.
 |          * 'all' : If all values are NA, drop that row or column.
 |      thresh : int, optional
 |          Require that many non-NA values.
 |      subset : array-like, optional
 |          Labels along other axis to consider, e.g. if you are dropping rows
 |          these would be a list of columns to include.
 |      inplace : bool, default False
 |          If True, do operation inplace and return None.
 |      
 |      Returns
 |      -------
 |      DataFrame
 |          DataFrame with NA entries dropped from it.
 |      
 |      See Also
 |      --------
 |      DataFrame.isna: Indicate missing values.
 |      DataFrame.notna : Indicate existing (non-missing) values.
 |      DataFrame.fillna : Replace missing values.
 |      Series.dropna : Drop missing values.
 |      Index.dropna : Drop missing indices.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({"name": ['Alfred', 'Batman', 'Catwoman'],
 |      ...                    "toy": [np.nan, 'Batmobile', 'Bullwhip'],
 |      ...                    "born": [pd.NaT, pd.Timestamp("1940-04-25"),
 |      ...                             pd.NaT]})
 |      >>> df
 |             name        toy       born
 |      0    Alfred        NaN        NaT
 |      1    Batman  Batmobile 1940-04-25
 |      2  Catwoman   Bullwhip        NaT
 |      
 |      Drop the rows where at least one element is missing.
 |      
 |      >>> df.dropna()
 |           name        toy       born
 |      1  Batman  Batmobile 1940-04-25
 |      
 |      Drop the columns where at least one element is missing.
 |      
 |      >>> df.dropna(axis='columns')
 |             name
 |      0    Alfred
 |      1    Batman
 |      2  Catwoman
 |      
 |      Drop the rows where all elements are missing.
 |      
 |      >>> df.dropna(how='all')
 |             name        toy       born
 |      0    Alfred        NaN        NaT
 |      1    Batman  Batmobile 1940-04-25
 |      2  Catwoman   Bullwhip        NaT
 |      
 |      Keep only the rows with at least 2 non-NA values.
 |      
 |      >>> df.dropna(thresh=2)
 |             name        toy       born
 |      1    Batman  Batmobile 1940-04-25
 |      2  Catwoman   Bullwhip        NaT
 |      
 |      Define in which columns to look for missing values.
 |      
 |      >>> df.dropna(subset=['name', 'born'])
 |             name        toy       born
 |      1    Batman  Batmobile 1940-04-25
 |      
 |      Keep the DataFrame with valid entries in the same variable.
 |      
 |      >>> df.dropna(inplace=True)
 |      >>> df
 |           name        toy       born
 |      1  Batman  Batmobile 1940-04-25
 |  
 |  duplicated(self, subset=None, keep='first')
 |      Return boolean Series denoting duplicate rows, optionally only
 |      considering certain columns
 |      
 |      Parameters
 |      ----------
 |      subset : column label or sequence of labels, optional
 |          Only consider certain columns for identifying duplicates, by
 |          default use all of the columns
 |      keep : {'first', 'last', False}, default 'first'
 |          - ``first`` : Mark duplicates as ``True`` except for the
 |            first occurrence.
 |          - ``last`` : Mark duplicates as ``True`` except for the
 |            last occurrence.
 |          - False : Mark all duplicates as ``True``.
 |      
 |      Returns
 |      -------
 |      duplicated : Series
 |  
 |  eq(self, other, axis='columns', level=None)
 |      Wrapper for flexible comparison methods eq
 |  
 |  eval(self, expr, inplace=False, **kwargs)
 |      Evaluate a string describing operations on DataFrame columns.
 |      
 |      Operates on columns only, not specific rows or elements.  This allows
 |      `eval` to run arbitrary code, which can make you vulnerable to code
 |      injection if you pass user input to this function.
 |      
 |      Parameters
 |      ----------
 |      expr : str
 |          The expression string to evaluate.
 |      inplace : bool, default False
 |          If the expression contains an assignment, whether to perform the
 |          operation inplace and mutate the existing DataFrame. Otherwise,
 |          a new DataFrame is returned.
 |      
 |          .. versionadded:: 0.18.0.
 |      kwargs : dict
 |          See the documentation for :func:`~pandas.eval` for complete details
 |          on the keyword arguments accepted by
 |          :meth:`~pandas.DataFrame.query`.
 |      
 |      Returns
 |      -------
 |      ndarray, scalar, or pandas object
 |          The result of the evaluation.
 |      
 |      See Also
 |      --------
 |      DataFrame.query : Evaluates a boolean expression to query the columns
 |          of a frame.
 |      DataFrame.assign : Can evaluate an expression or function to create new
 |          values for a column.
 |      pandas.eval : Evaluate a Python expression as a string using various
 |          backends.
 |      
 |      Notes
 |      -----
 |      For more details see the API documentation for :func:`~pandas.eval`.
 |      For detailed examples see :ref:`enhancing performance with eval
 |      <enhancingperf.eval>`.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'A': range(1, 6), 'B': range(10, 0, -2)})
 |      >>> df
 |         A   B
 |      0  1  10
 |      1  2   8
 |      2  3   6
 |      3  4   4
 |      4  5   2
 |      >>> df.eval('A + B')
 |      0    11
 |      1    10
 |      2     9
 |      3     8
 |      4     7
 |      dtype: int64
 |      
 |      Assignment is allowed though by default the original DataFrame is not
 |      modified.
 |      
 |      >>> df.eval('C = A + B')
 |         A   B   C
 |      0  1  10  11
 |      1  2   8  10
 |      2  3   6   9
 |      3  4   4   8
 |      4  5   2   7
 |      >>> df
 |         A   B
 |      0  1  10
 |      1  2   8
 |      2  3   6
 |      3  4   4
 |      4  5   2
 |      
 |      Use ``inplace=True`` to modify the original DataFrame.
 |      
 |      >>> df.eval('C = A + B', inplace=True)
 |      >>> df
 |         A   B   C
 |      0  1  10  11
 |      1  2   8  10
 |      2  3   6   9
 |      3  4   4   8
 |      4  5   2   7
 |  
 |  ewm(self, com=None, span=None, halflife=None, alpha=None, min_periods=0, adjust=True, ignore_na=False, axis=0)
 |      Provides exponential weighted functions
 |      
 |      .. versionadded:: 0.18.0
 |      
 |      Parameters
 |      ----------
 |      com : float, optional
 |          Specify decay in terms of center of mass,
 |          :math:`\alpha = 1 / (1 + com),\text{ for } com \geq 0`
 |      span : float, optional
 |          Specify decay in terms of span,
 |          :math:`\alpha = 2 / (span + 1),\text{ for } span \geq 1`
 |      halflife : float, optional
 |          Specify decay in terms of half-life,
 |          :math:`\alpha = 1 - exp(log(0.5) / halflife),\text{ for } halflife > 0`
 |      alpha : float, optional
 |          Specify smoothing factor :math:`\alpha` directly,
 |          :math:`0 < \alpha \leq 1`
 |      
 |          .. versionadded:: 0.18.0
 |      
 |      min_periods : int, default 0
 |          Minimum number of observations in window required to have a value
 |          (otherwise result is NA).
 |      adjust : boolean, default True
 |          Divide by decaying adjustment factor in beginning periods to account
 |          for imbalance in relative weightings (viewing EWMA as a moving average)
 |      ignore_na : boolean, default False
 |          Ignore missing values when calculating weights;
 |          specify True to reproduce pre-0.15.0 behavior
 |      
 |      Returns
 |      -------
 |      a Window sub-classed for the particular operation
 |      
 |      Examples
 |      --------
 |      
 |      >>> df = DataFrame({'B': [0, 1, 2, np.nan, 4]})
 |           B
 |      0  0.0
 |      1  1.0
 |      2  2.0
 |      3  NaN
 |      4  4.0
 |      
 |      >>> df.ewm(com=0.5).mean()
 |                B
 |      0  0.000000
 |      1  0.750000
 |      2  1.615385
 |      3  1.615385
 |      4  3.670213
 |      
 |      Notes
 |      -----
 |      Exactly one of center of mass, span, half-life, and alpha must be provided.
 |      Allowed values and relationship between the parameters are specified in the
 |      parameter descriptions above; see the link at the end of this section for
 |      a detailed explanation.
 |      
 |      When adjust is True (default), weighted averages are calculated using
 |      weights (1-alpha)**(n-1), (1-alpha)**(n-2), ..., 1-alpha, 1.
 |      
 |      When adjust is False, weighted averages are calculated recursively as:
 |         weighted_average[0] = arg[0];
 |         weighted_average[i] = (1-alpha)*weighted_average[i-1] + alpha*arg[i].
 |      
 |      When ignore_na is False (default), weights are based on absolute positions.
 |      For example, the weights of x and y used in calculating the final weighted
 |      average of [x, None, y] are (1-alpha)**2 and 1 (if adjust is True), and
 |      (1-alpha)**2 and alpha (if adjust is False).
 |      
 |      When ignore_na is True (reproducing pre-0.15.0 behavior), weights are based
 |      on relative positions. For example, the weights of x and y used in
 |      calculating the final weighted average of [x, None, y] are 1-alpha and 1
 |      (if adjust is True), and 1-alpha and alpha (if adjust is False).
 |      
 |      More details can be found at
 |      http://pandas.pydata.org/pandas-docs/stable/computation.html#exponentially-weighted-windows
 |      
 |      See Also
 |      --------
 |      rolling : Provides rolling window calculations
 |      expanding : Provides expanding transformations.
 |  
 |  expanding(self, min_periods=1, center=False, axis=0)
 |      Provides expanding transformations.
 |      
 |      .. versionadded:: 0.18.0
 |      
 |      Parameters
 |      ----------
 |      min_periods : int, default 1
 |          Minimum number of observations in window required to have a value
 |          (otherwise result is NA).
 |      center : boolean, default False
 |          Set the labels at the center of the window.
 |      axis : int or string, default 0
 |      
 |      Returns
 |      -------
 |      a Window sub-classed for the particular operation
 |      
 |      Examples
 |      --------
 |      
 |      >>> df = DataFrame({'B': [0, 1, 2, np.nan, 4]})
 |           B
 |      0  0.0
 |      1  1.0
 |      2  2.0
 |      3  NaN
 |      4  4.0
 |      
 |      >>> df.expanding(2).sum()
 |           B
 |      0  NaN
 |      1  1.0
 |      2  3.0
 |      3  3.0
 |      4  7.0
 |      
 |      Notes
 |      -----
 |      By default, the result is set to the right edge of the window. This can be
 |      changed to the center of the window by setting ``center=True``.
 |      
 |      See Also
 |      --------
 |      rolling : Provides rolling window calculations
 |      ewm : Provides exponential weighted functions
 |  
 |  fillna(self, value=None, method=None, axis=None, inplace=False, limit=None, downcast=None, **kwargs)
 |      Fill NA/NaN values using the specified method
 |      
 |      Parameters
 |      ----------
 |      value : scalar, dict, Series, or DataFrame
 |          Value to use to fill holes (e.g. 0), alternately a
 |          dict/Series/DataFrame of values specifying which value to use for
 |          each index (for a Series) or column (for a DataFrame). (values not
 |          in the dict/Series/DataFrame will not be filled). This value cannot
 |          be a list.
 |      method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
 |          Method to use for filling holes in reindexed Series
 |          pad / ffill: propagate last valid observation forward to next valid
 |          backfill / bfill: use NEXT valid observation to fill gap
 |      axis : {0 or 'index', 1 or 'columns'}
 |      inplace : boolean, default False
 |          If True, fill in place. Note: this will modify any
 |          other views on this object, (e.g. a no-copy slice for a column in a
 |          DataFrame).
 |      limit : int, default None
 |          If method is specified, this is the maximum number of consecutive
 |          NaN values to forward/backward fill. In other words, if there is
 |          a gap with more than this number of consecutive NaNs, it will only
 |          be partially filled. If method is not specified, this is the
 |          maximum number of entries along the entire axis where NaNs will be
 |          filled. Must be greater than 0 if not None.
 |      downcast : dict, default is None
 |          a dict of item->dtype of what to downcast if possible,
 |          or the string 'infer' which will try to downcast to an appropriate
 |          equal type (e.g. float64 to int64 if possible)
 |      
 |      See Also
 |      --------
 |      interpolate : Fill NaN values using interpolation.
 |      reindex, asfreq
 |      
 |      Returns
 |      -------
 |      filled : DataFrame
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame([[np.nan, 2, np.nan, 0],
 |      ...                    [3, 4, np.nan, 1],
 |      ...                    [np.nan, np.nan, np.nan, 5],
 |      ...                    [np.nan, 3, np.nan, 4]],
 |      ...                    columns=list('ABCD'))
 |      >>> df
 |           A    B   C  D
 |      0  NaN  2.0 NaN  0
 |      1  3.0  4.0 NaN  1
 |      2  NaN  NaN NaN  5
 |      3  NaN  3.0 NaN  4
 |      
 |      Replace all NaN elements with 0s.
 |      
 |      >>> df.fillna(0)
 |          A   B   C   D
 |      0   0.0 2.0 0.0 0
 |      1   3.0 4.0 0.0 1
 |      2   0.0 0.0 0.0 5
 |      3   0.0 3.0 0.0 4
 |      
 |      We can also propagate non-null values forward or backward.
 |      
 |      >>> df.fillna(method='ffill')
 |          A   B   C   D
 |      0   NaN 2.0 NaN 0
 |      1   3.0 4.0 NaN 1
 |      2   3.0 4.0 NaN 5
 |      3   3.0 3.0 NaN 4
 |      
 |      Replace all NaN elements in column 'A', 'B', 'C', and 'D', with 0, 1,
 |      2, and 3 respectively.
 |      
 |      >>> values = {'A': 0, 'B': 1, 'C': 2, 'D': 3}
 |      >>> df.fillna(value=values)
 |          A   B   C   D
 |      0   0.0 2.0 2.0 0
 |      1   3.0 4.0 2.0 1
 |      2   0.0 1.0 2.0 5
 |      3   0.0 3.0 2.0 4
 |      
 |      Only replace the first NaN element.
 |      
 |      >>> df.fillna(value=values, limit=1)
 |          A   B   C   D
 |      0   0.0 2.0 2.0 0
 |      1   3.0 4.0 NaN 1
 |      2   NaN 1.0 NaN 5
 |      3   NaN 3.0 NaN 4
 |  
 |  floordiv(self, other, axis='columns', level=None, fill_value=None)
 |      Integer division of dataframe and other, element-wise (binary operator `floordiv`).
 |      
 |      Equivalent to ``dataframe // other``, but with support to substitute a fill_value for
 |      missing data in one of the inputs.
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      None
 |      
 |      See also
 |      --------
 |      DataFrame.rfloordiv
 |  
 |  ge(self, other, axis='columns', level=None)
 |      Wrapper for flexible comparison methods ge
 |  
 |  get_value(self, index, col, takeable=False)
 |      Quickly retrieve single value at passed column and index
 |      
 |      .. deprecated:: 0.21.0
 |          Use .at[] or .iat[] accessors instead.
 |      
 |      Parameters
 |      ----------
 |      index : row label
 |      col : column label
 |      takeable : interpret the index/col as indexers, default False
 |      
 |      Returns
 |      -------
 |      value : scalar value
 |  
 |  gt(self, other, axis='columns', level=None)
 |      Wrapper for flexible comparison methods gt
 |  
 |  hist = hist_frame(data, column=None, by=None, grid=True, xlabelsize=None, xrot=None, ylabelsize=None, yrot=None, ax=None, sharex=False, sharey=False, figsize=None, layout=None, bins=10, **kwds)
 |      Make a histogram of the DataFrame's.
 |      
 |      A `histogram`_ is a representation of the distribution of data.
 |      This function calls :meth:`matplotlib.pyplot.hist`, on each series in
 |      the DataFrame, resulting in one histogram per column.
 |      
 |      .. _histogram: https://en.wikipedia.org/wiki/Histogram
 |      
 |      Parameters
 |      ----------
 |      data : DataFrame
 |          The pandas object holding the data.
 |      column : string or sequence
 |          If passed, will be used to limit data to a subset of columns.
 |      by : object, optional
 |          If passed, then used to form histograms for separate groups.
 |      grid : boolean, default True
 |          Whether to show axis grid lines.
 |      xlabelsize : int, default None
 |          If specified changes the x-axis label size.
 |      xrot : float, default None
 |          Rotation of x axis labels. For example, a value of 90 displays the
 |          x labels rotated 90 degrees clockwise.
 |      ylabelsize : int, default None
 |          If specified changes the y-axis label size.
 |      yrot : float, default None
 |          Rotation of y axis labels. For example, a value of 90 displays the
 |          y labels rotated 90 degrees clockwise.
 |      ax : Matplotlib axes object, default None
 |          The axes to plot the histogram on.
 |      sharex : boolean, default True if ax is None else False
 |          In case subplots=True, share x axis and set some x axis labels to
 |          invisible; defaults to True if ax is None otherwise False if an ax
 |          is passed in.
 |          Note that passing in both an ax and sharex=True will alter all x axis
 |          labels for all subplots in a figure.
 |      sharey : boolean, default False
 |          In case subplots=True, share y axis and set some y axis labels to
 |          invisible.
 |      figsize : tuple
 |          The size in inches of the figure to create. Uses the value in
 |          `matplotlib.rcParams` by default.
 |      layout : tuple, optional
 |          Tuple of (rows, columns) for the layout of the histograms.
 |      bins : integer or sequence, default 10
 |          Number of histogram bins to be used. If an integer is given, bins + 1
 |          bin edges are calculated and returned. If bins is a sequence, gives
 |          bin edges, including left edge of first bin and right edge of last
 |          bin. In this case, bins is returned unmodified.
 |      **kwds
 |          All other plotting keyword arguments to be passed to
 |          :meth:`matplotlib.pyplot.hist`.
 |      
 |      Returns
 |      -------
 |      axes : matplotlib.AxesSubplot or numpy.ndarray of them
 |      
 |      See Also
 |      --------
 |      matplotlib.pyplot.hist : Plot a histogram using matplotlib.
 |      
 |      Examples
 |      --------
 |      
 |      .. plot::
 |          :context: close-figs
 |      
 |          This example draws a histogram based on the length and width of
 |          some animals, displayed in three bins
 |      
 |          >>> df = pd.DataFrame({
 |          ...     'length': [1.5, 0.5, 1.2, 0.9, 3],
 |          ...     'width': [0.7, 0.2, 0.15, 0.2, 1.1]
 |          ...     }, index= ['pig', 'rabbit', 'duck', 'chicken', 'horse'])
 |          >>> hist = df.hist(bins=3)
 |  
 |  idxmax(self, axis=0, skipna=True)
 |      Return index of first occurrence of maximum over requested axis.
 |      NA/null values are excluded.
 |      
 |      Parameters
 |      ----------
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          0 or 'index' for row-wise, 1 or 'columns' for column-wise
 |      skipna : boolean, default True
 |          Exclude NA/null values. If an entire row/column is NA, the result
 |          will be NA.
 |      
 |      Raises
 |      ------
 |      ValueError
 |          * If the row/column is empty
 |      
 |      Returns
 |      -------
 |      idxmax : Series
 |      
 |      Notes
 |      -----
 |      This method is the DataFrame version of ``ndarray.argmax``.
 |      
 |      See Also
 |      --------
 |      Series.idxmax
 |  
 |  idxmin(self, axis=0, skipna=True)
 |      Return index of first occurrence of minimum over requested axis.
 |      NA/null values are excluded.
 |      
 |      Parameters
 |      ----------
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          0 or 'index' for row-wise, 1 or 'columns' for column-wise
 |      skipna : boolean, default True
 |          Exclude NA/null values. If an entire row/column is NA, the result
 |          will be NA.
 |      
 |      Raises
 |      ------
 |      ValueError
 |          * If the row/column is empty
 |      
 |      Returns
 |      -------
 |      idxmin : Series
 |      
 |      Notes
 |      -----
 |      This method is the DataFrame version of ``ndarray.argmin``.
 |      
 |      See Also
 |      --------
 |      Series.idxmin
 |  
 |  info(self, verbose=None, buf=None, max_cols=None, memory_usage=None, null_counts=None)
 |      Print a concise summary of a DataFrame.
 |      
 |      This method prints information about a DataFrame including
 |      the index dtype and column dtypes, non-null values and memory usage.
 |      
 |      Parameters
 |      ----------
 |      verbose : bool, optional
 |          Whether to print the full summary. By default, the setting in
 |          ``pandas.options.display.max_info_columns`` is followed.
 |      buf : writable buffer, defaults to sys.stdout
 |          Where to send the output. By default, the output is printed to
 |          sys.stdout. Pass a writable buffer if you need to further process
 |          the output.
 |      max_cols : int, optional
 |          When to switch from the verbose to the truncated output. If the
 |          DataFrame has more than `max_cols` columns, the truncated output
 |          is used. By default, the setting in
 |          ``pandas.options.display.max_info_columns`` is used.
 |      memory_usage : bool, str, optional
 |          Specifies whether total memory usage of the DataFrame
 |          elements (including the index) should be displayed. By default,
 |          this follows the ``pandas.options.display.memory_usage`` setting.
 |      
 |          True always show memory usage. False never shows memory usage.
 |          A value of 'deep' is equivalent to "True with deep introspection".
 |          Memory usage is shown in human-readable units (base-2
 |          representation). Without deep introspection a memory estimation is
 |          made based in column dtype and number of rows assuming values
 |          consume the same memory amount for corresponding dtypes. With deep
 |          memory introspection, a real memory usage calculation is performed
 |          at the cost of computational resources.
 |      null_counts : bool, optional
 |          Whether to show the non-null counts. By default, this is shown
 |          only if the frame is smaller than
 |          ``pandas.options.display.max_info_rows`` and
 |          ``pandas.options.display.max_info_columns``. A value of True always
 |          shows the counts, and False never shows the counts.
 |      
 |      Returns
 |      -------
 |      None
 |          This method prints a summary of a DataFrame and returns None.
 |      
 |      See Also
 |      --------
 |      DataFrame.describe: Generate descriptive statistics of DataFrame
 |          columns.
 |      DataFrame.memory_usage: Memory usage of DataFrame columns.
 |      
 |      Examples
 |      --------
 |      >>> int_values = [1, 2, 3, 4, 5]
 |      >>> text_values = ['alpha', 'beta', 'gamma', 'delta', 'epsilon']
 |      >>> float_values = [0.0, 0.25, 0.5, 0.75, 1.0]
 |      >>> df = pd.DataFrame({"int_col": int_values, "text_col": text_values,
 |      ...                   "float_col": float_values})
 |      >>> df
 |         int_col text_col  float_col
 |      0        1    alpha       0.00
 |      1        2     beta       0.25
 |      2        3    gamma       0.50
 |      3        4    delta       0.75
 |      4        5  epsilon       1.00
 |      
 |      Prints information of all columns:
 |      
 |      >>> df.info(verbose=True)
 |      <class 'pandas.core.frame.DataFrame'>
 |      RangeIndex: 5 entries, 0 to 4
 |      Data columns (total 3 columns):
 |      int_col      5 non-null int64
 |      text_col     5 non-null object
 |      float_col    5 non-null float64
 |      dtypes: float64(1), int64(1), object(1)
 |      memory usage: 200.0+ bytes
 |      
 |      Prints a summary of columns count and its dtypes but not per column
 |      information:
 |      
 |      >>> df.info(verbose=False)
 |      <class 'pandas.core.frame.DataFrame'>
 |      RangeIndex: 5 entries, 0 to 4
 |      Columns: 3 entries, int_col to float_col
 |      dtypes: float64(1), int64(1), object(1)
 |      memory usage: 200.0+ bytes
 |      
 |      Pipe output of DataFrame.info to buffer instead of sys.stdout, get
 |      buffer content and writes to a text file:
 |      
 |      >>> import io
 |      >>> buffer = io.StringIO()
 |      >>> df.info(buf=buffer)
 |      >>> s = buffer.getvalue()
 |      >>> with open("df_info.txt", "w", encoding="utf-8") as f:
 |      ...     f.write(s)
 |      260
 |      
 |      The `memory_usage` parameter allows deep introspection mode, specially
 |      useful for big DataFrames and fine-tune memory optimization:
 |      
 |      >>> random_strings_array = np.random.choice(['a', 'b', 'c'], 10 ** 6)
 |      >>> df = pd.DataFrame({
 |      ...     'column_1': np.random.choice(['a', 'b', 'c'], 10 ** 6),
 |      ...     'column_2': np.random.choice(['a', 'b', 'c'], 10 ** 6),
 |      ...     'column_3': np.random.choice(['a', 'b', 'c'], 10 ** 6)
 |      ... })
 |      >>> df.info()
 |      <class 'pandas.core.frame.DataFrame'>
 |      RangeIndex: 1000000 entries, 0 to 999999
 |      Data columns (total 3 columns):
 |      column_1    1000000 non-null object
 |      column_2    1000000 non-null object
 |      column_3    1000000 non-null object
 |      dtypes: object(3)
 |      memory usage: 22.9+ MB
 |      
 |      >>> df.info(memory_usage='deep')
 |      <class 'pandas.core.frame.DataFrame'>
 |      RangeIndex: 1000000 entries, 0 to 999999
 |      Data columns (total 3 columns):
 |      column_1    1000000 non-null object
 |      column_2    1000000 non-null object
 |      column_3    1000000 non-null object
 |      dtypes: object(3)
 |      memory usage: 188.8 MB
 |  
 |  insert(self, loc, column, value, allow_duplicates=False)
 |      Insert column into DataFrame at specified location.
 |      
 |      Raises a ValueError if `column` is already contained in the DataFrame,
 |      unless `allow_duplicates` is set to True.
 |      
 |      Parameters
 |      ----------
 |      loc : int
 |          Insertion index. Must verify 0 <= loc <= len(columns)
 |      column : string, number, or hashable object
 |          label of the inserted column
 |      value : int, Series, or array-like
 |      allow_duplicates : bool, optional
 |  
 |  isin(self, values)
 |      Return boolean DataFrame showing whether each element in the
 |      DataFrame is contained in values.
 |      
 |      Parameters
 |      ----------
 |      values : iterable, Series, DataFrame or dictionary
 |          The result will only be true at a location if all the
 |          labels match. If `values` is a Series, that's the index. If
 |          `values` is a dictionary, the keys must be the column names,
 |          which must match. If `values` is a DataFrame,
 |          then both the index and column labels must match.
 |      
 |      Returns
 |      -------
 |      
 |      DataFrame of booleans
 |      
 |      Examples
 |      --------
 |      When ``values`` is a list:
 |      
 |      >>> df = pd.DataFrame({'A': [1, 2, 3], 'B': ['a', 'b', 'f']})
 |      >>> df.isin([1, 3, 12, 'a'])
 |             A      B
 |      0   True   True
 |      1  False  False
 |      2   True  False
 |      
 |      When ``values`` is a dict:
 |      
 |      >>> df = pd.DataFrame({'A': [1, 2, 3], 'B': [1, 4, 7]})
 |      >>> df.isin({'A': [1, 3], 'B': [4, 7, 12]})
 |             A      B
 |      0   True  False  # Note that B didn't match the 1 here.
 |      1  False   True
 |      2   True   True
 |      
 |      When ``values`` is a Series or DataFrame:
 |      
 |      >>> df = pd.DataFrame({'A': [1, 2, 3], 'B': ['a', 'b', 'f']})
 |      >>> other = DataFrame({'A': [1, 3, 3, 2], 'B': ['e', 'f', 'f', 'e']})
 |      >>> df.isin(other)
 |             A      B
 |      0   True  False
 |      1  False  False  # Column A in `other` has a 3, but not at index 1.
 |      2   True   True
 |  
 |  isna(self)
 |      Detect missing values.
 |      
 |      Return a boolean same-sized object indicating if the values are NA.
 |      NA values, such as None or :attr:`numpy.NaN`, gets mapped to True
 |      values.
 |      Everything else gets mapped to False values. Characters such as empty
 |      strings ``''`` or :attr:`numpy.inf` are not considered NA values
 |      (unless you set ``pandas.options.mode.use_inf_as_na = True``).
 |      
 |      Returns
 |      -------
 |      DataFrame
 |          Mask of bool values for each element in DataFrame that
 |          indicates whether an element is not an NA value.
 |      
 |      See Also
 |      --------
 |      DataFrame.isnull : alias of isna
 |      DataFrame.notna : boolean inverse of isna
 |      DataFrame.dropna : omit axes labels with missing values
 |      isna : top-level isna
 |      
 |      Examples
 |      --------
 |      Show which entries in a DataFrame are NA.
 |      
 |      >>> df = pd.DataFrame({'age': [5, 6, np.NaN],
 |      ...                    'born': [pd.NaT, pd.Timestamp('1939-05-27'),
 |      ...                             pd.Timestamp('1940-04-25')],
 |      ...                    'name': ['Alfred', 'Batman', ''],
 |      ...                    'toy': [None, 'Batmobile', 'Joker']})
 |      >>> df
 |         age       born    name        toy
 |      0  5.0        NaT  Alfred       None
 |      1  6.0 1939-05-27  Batman  Batmobile
 |      2  NaN 1940-04-25              Joker
 |      
 |      >>> df.isna()
 |           age   born   name    toy
 |      0  False   True  False   True
 |      1  False  False  False  False
 |      2   True  False  False  False
 |      
 |      Show which entries in a Series are NA.
 |      
 |      >>> ser = pd.Series([5, 6, np.NaN])
 |      >>> ser
 |      0    5.0
 |      1    6.0
 |      2    NaN
 |      dtype: float64
 |      
 |      >>> ser.isna()
 |      0    False
 |      1    False
 |      2     True
 |      dtype: bool
 |  
 |  isnull(self)
 |      Detect missing values.
 |      
 |      Return a boolean same-sized object indicating if the values are NA.
 |      NA values, such as None or :attr:`numpy.NaN`, gets mapped to True
 |      values.
 |      Everything else gets mapped to False values. Characters such as empty
 |      strings ``''`` or :attr:`numpy.inf` are not considered NA values
 |      (unless you set ``pandas.options.mode.use_inf_as_na = True``).
 |      
 |      Returns
 |      -------
 |      DataFrame
 |          Mask of bool values for each element in DataFrame that
 |          indicates whether an element is not an NA value.
 |      
 |      See Also
 |      --------
 |      DataFrame.isnull : alias of isna
 |      DataFrame.notna : boolean inverse of isna
 |      DataFrame.dropna : omit axes labels with missing values
 |      isna : top-level isna
 |      
 |      Examples
 |      --------
 |      Show which entries in a DataFrame are NA.
 |      
 |      >>> df = pd.DataFrame({'age': [5, 6, np.NaN],
 |      ...                    'born': [pd.NaT, pd.Timestamp('1939-05-27'),
 |      ...                             pd.Timestamp('1940-04-25')],
 |      ...                    'name': ['Alfred', 'Batman', ''],
 |      ...                    'toy': [None, 'Batmobile', 'Joker']})
 |      >>> df
 |         age       born    name        toy
 |      0  5.0        NaT  Alfred       None
 |      1  6.0 1939-05-27  Batman  Batmobile
 |      2  NaN 1940-04-25              Joker
 |      
 |      >>> df.isna()
 |           age   born   name    toy
 |      0  False   True  False   True
 |      1  False  False  False  False
 |      2   True  False  False  False
 |      
 |      Show which entries in a Series are NA.
 |      
 |      >>> ser = pd.Series([5, 6, np.NaN])
 |      >>> ser
 |      0    5.0
 |      1    6.0
 |      2    NaN
 |      dtype: float64
 |      
 |      >>> ser.isna()
 |      0    False
 |      1    False
 |      2     True
 |      dtype: bool
 |  
 |  items = iteritems(self)
 |  
 |  iteritems(self)
 |      Iterator over (column name, Series) pairs.
 |      
 |      See also
 |      --------
 |      iterrows : Iterate over DataFrame rows as (index, Series) pairs.
 |      itertuples : Iterate over DataFrame rows as namedtuples of the values.
 |  
 |  iterrows(self)
 |      Iterate over DataFrame rows as (index, Series) pairs.
 |      
 |      Notes
 |      -----
 |      
 |      1. Because ``iterrows`` returns a Series for each row,
 |         it does **not** preserve dtypes across the rows (dtypes are
 |         preserved across columns for DataFrames). For example,
 |      
 |         >>> df = pd.DataFrame([[1, 1.5]], columns=['int', 'float'])
 |         >>> row = next(df.iterrows())[1]
 |         >>> row
 |         int      1.0
 |         float    1.5
 |         Name: 0, dtype: float64
 |         >>> print(row['int'].dtype)
 |         float64
 |         >>> print(df['int'].dtype)
 |         int64
 |      
 |         To preserve dtypes while iterating over the rows, it is better
 |         to use :meth:`itertuples` which returns namedtuples of the values
 |         and which is generally faster than ``iterrows``.
 |      
 |      2. You should **never modify** something you are iterating over.
 |         This is not guaranteed to work in all cases. Depending on the
 |         data types, the iterator returns a copy and not a view, and writing
 |         to it will have no effect.
 |      
 |      Returns
 |      -------
 |      it : generator
 |          A generator that iterates over the rows of the frame.
 |      
 |      See also
 |      --------
 |      itertuples : Iterate over DataFrame rows as namedtuples of the values.
 |      iteritems : Iterate over (column name, Series) pairs.
 |  
 |  itertuples(self, index=True, name='Pandas')
 |      Iterate over DataFrame rows as namedtuples, with index value as first
 |      element of the tuple.
 |      
 |      Parameters
 |      ----------
 |      index : boolean, default True
 |          If True, return the index as the first element of the tuple.
 |      name : string, default "Pandas"
 |          The name of the returned namedtuples or None to return regular
 |          tuples.
 |      
 |      Notes
 |      -----
 |      The column names will be renamed to positional names if they are
 |      invalid Python identifiers, repeated, or start with an underscore.
 |      With a large number of columns (>255), regular tuples are returned.
 |      
 |      See also
 |      --------
 |      iterrows : Iterate over DataFrame rows as (index, Series) pairs.
 |      iteritems : Iterate over (column name, Series) pairs.
 |      
 |      Examples
 |      --------
 |      
 |      >>> df = pd.DataFrame({'col1': [1, 2], 'col2': [0.1, 0.2]},
 |                            index=['a', 'b'])
 |      >>> df
 |         col1  col2
 |      a     1   0.1
 |      b     2   0.2
 |      >>> for row in df.itertuples():
 |      ...     print(row)
 |      ...
 |      Pandas(Index='a', col1=1, col2=0.10000000000000001)
 |      Pandas(Index='b', col1=2, col2=0.20000000000000001)
 |  
 |  join(self, other, on=None, how='left', lsuffix='', rsuffix='', sort=False)
 |      Join columns with other DataFrame either on index or on a key
 |      column. Efficiently Join multiple DataFrame objects by index at once by
 |      passing a list.
 |      
 |      Parameters
 |      ----------
 |      other : DataFrame, Series with name field set, or list of DataFrame
 |          Index should be similar to one of the columns in this one. If a
 |          Series is passed, its name attribute must be set, and that will be
 |          used as the column name in the resulting joined DataFrame
 |      on : name, tuple/list of names, or array-like
 |          Column or index level name(s) in the caller to join on the index
 |          in `other`, otherwise joins index-on-index. If multiple
 |          values given, the `other` DataFrame must have a MultiIndex. Can
 |          pass an array as the join key if it is not already contained in
 |          the calling DataFrame. Like an Excel VLOOKUP operation
 |      how : {'left', 'right', 'outer', 'inner'}, default: 'left'
 |          How to handle the operation of the two objects.
 |      
 |          * left: use calling frame's index (or column if on is specified)
 |          * right: use other frame's index
 |          * outer: form union of calling frame's index (or column if on is
 |            specified) with other frame's index, and sort it
 |            lexicographically
 |          * inner: form intersection of calling frame's index (or column if
 |            on is specified) with other frame's index, preserving the order
 |            of the calling's one
 |      lsuffix : string
 |          Suffix to use from left frame's overlapping columns
 |      rsuffix : string
 |          Suffix to use from right frame's overlapping columns
 |      sort : boolean, default False
 |          Order result DataFrame lexicographically by the join key. If False,
 |          the order of the join key depends on the join type (how keyword)
 |      
 |      Notes
 |      -----
 |      on, lsuffix, and rsuffix options are not supported when passing a list
 |      of DataFrame objects
 |      
 |      Support for specifying index levels as the `on` parameter was added
 |      in version 0.23.0
 |      
 |      Examples
 |      --------
 |      >>> caller = pd.DataFrame({'key': ['K0', 'K1', 'K2', 'K3', 'K4', 'K5'],
 |      ...                        'A': ['A0', 'A1', 'A2', 'A3', 'A4', 'A5']})
 |      
 |      >>> caller
 |          A key
 |      0  A0  K0
 |      1  A1  K1
 |      2  A2  K2
 |      3  A3  K3
 |      4  A4  K4
 |      5  A5  K5
 |      
 |      >>> other = pd.DataFrame({'key': ['K0', 'K1', 'K2'],
 |      ...                       'B': ['B0', 'B1', 'B2']})
 |      
 |      >>> other
 |          B key
 |      0  B0  K0
 |      1  B1  K1
 |      2  B2  K2
 |      
 |      Join DataFrames using their indexes.
 |      
 |      >>> caller.join(other, lsuffix='_caller', rsuffix='_other')
 |      
 |      >>>     A key_caller    B key_other
 |          0  A0         K0   B0        K0
 |          1  A1         K1   B1        K1
 |          2  A2         K2   B2        K2
 |          3  A3         K3  NaN       NaN
 |          4  A4         K4  NaN       NaN
 |          5  A5         K5  NaN       NaN
 |      
 |      
 |      If we want to join using the key columns, we need to set key to be
 |      the index in both caller and other. The joined DataFrame will have
 |      key as its index.
 |      
 |      >>> caller.set_index('key').join(other.set_index('key'))
 |      
 |      >>>      A    B
 |          key
 |          K0   A0   B0
 |          K1   A1   B1
 |          K2   A2   B2
 |          K3   A3  NaN
 |          K4   A4  NaN
 |          K5   A5  NaN
 |      
 |      Another option to join using the key columns is to use the on
 |      parameter. DataFrame.join always uses other's index but we can use any
 |      column in the caller. This method preserves the original caller's
 |      index in the result.
 |      
 |      >>> caller.join(other.set_index('key'), on='key')
 |      
 |      >>>     A key    B
 |          0  A0  K0   B0
 |          1  A1  K1   B1
 |          2  A2  K2   B2
 |          3  A3  K3  NaN
 |          4  A4  K4  NaN
 |          5  A5  K5  NaN
 |      
 |      
 |      See also
 |      --------
 |      DataFrame.merge : For column(s)-on-columns(s) operations
 |      
 |      Returns
 |      -------
 |      joined : DataFrame
 |  
 |  kurt(self, axis=None, skipna=None, level=None, numeric_only=None, **kwargs)
 |      Return unbiased kurtosis over requested axis using Fisher's definition of
 |      kurtosis (kurtosis of normal == 0.0). Normalized by N-1
 |      
 |      
 |      Parameters
 |      ----------
 |      axis : {index (0), columns (1)}
 |      skipna : boolean, default True
 |          Exclude NA/null values when computing the result.
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series
 |      numeric_only : boolean, default None
 |          Include only float, int, boolean columns. If None, will attempt to use
 |          everything, then use only numeric data. Not implemented for Series.
 |      
 |      Returns
 |      -------
 |      kurt : Series or DataFrame (if level specified)
 |  
 |  kurtosis = kurt(self, axis=None, skipna=None, level=None, numeric_only=None, **kwargs)
 |  
 |  le(self, other, axis='columns', level=None)
 |      Wrapper for flexible comparison methods le
 |  
 |  lookup(self, row_labels, col_labels)
 |      Label-based "fancy indexing" function for DataFrame.
 |      Given equal-length arrays of row and column labels, return an
 |      array of the values corresponding to each (row, col) pair.
 |      
 |      Parameters
 |      ----------
 |      row_labels : sequence
 |          The row labels to use for lookup
 |      col_labels : sequence
 |          The column labels to use for lookup
 |      
 |      Notes
 |      -----
 |      Akin to::
 |      
 |          result = []
 |          for row, col in zip(row_labels, col_labels):
 |              result.append(df.get_value(row, col))
 |      
 |      Examples
 |      --------
 |      values : ndarray
 |          The found values
 |  
 |  lt(self, other, axis='columns', level=None)
 |      Wrapper for flexible comparison methods lt
 |  
 |  mad(self, axis=None, skipna=None, level=None)
 |      Return the mean absolute deviation of the values for the requested axis
 |      
 |      Parameters
 |      ----------
 |      axis : {index (0), columns (1)}
 |      skipna : boolean, default True
 |          Exclude NA/null values when computing the result.
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series
 |      numeric_only : boolean, default None
 |          Include only float, int, boolean columns. If None, will attempt to use
 |          everything, then use only numeric data. Not implemented for Series.
 |      
 |      Returns
 |      -------
 |      mad : Series or DataFrame (if level specified)
 |  
 |  max(self, axis=None, skipna=None, level=None, numeric_only=None, **kwargs)
 |      This method returns the maximum of the values in the object.
 |                  If you want the *index* of the maximum, use ``idxmax``. This is
 |                  the equivalent of the ``numpy.ndarray`` method ``argmax``.
 |      
 |      Parameters
 |      ----------
 |      axis : {index (0), columns (1)}
 |      skipna : boolean, default True
 |          Exclude NA/null values when computing the result.
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series
 |      numeric_only : boolean, default None
 |          Include only float, int, boolean columns. If None, will attempt to use
 |          everything, then use only numeric data. Not implemented for Series.
 |      
 |      Returns
 |      -------
 |      max : Series or DataFrame (if level specified)
 |  
 |  mean(self, axis=None, skipna=None, level=None, numeric_only=None, **kwargs)
 |      Return the mean of the values for the requested axis
 |      
 |      Parameters
 |      ----------
 |      axis : {index (0), columns (1)}
 |      skipna : boolean, default True
 |          Exclude NA/null values when computing the result.
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series
 |      numeric_only : boolean, default None
 |          Include only float, int, boolean columns. If None, will attempt to use
 |          everything, then use only numeric data. Not implemented for Series.
 |      
 |      Returns
 |      -------
 |      mean : Series or DataFrame (if level specified)
 |  
 |  median(self, axis=None, skipna=None, level=None, numeric_only=None, **kwargs)
 |      Return the median of the values for the requested axis
 |      
 |      Parameters
 |      ----------
 |      axis : {index (0), columns (1)}
 |      skipna : boolean, default True
 |          Exclude NA/null values when computing the result.
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series
 |      numeric_only : boolean, default None
 |          Include only float, int, boolean columns. If None, will attempt to use
 |          everything, then use only numeric data. Not implemented for Series.
 |      
 |      Returns
 |      -------
 |      median : Series or DataFrame (if level specified)
 |  
 |  melt(self, id_vars=None, value_vars=None, var_name=None, value_name='value', col_level=None)
 |      "Unpivots" a DataFrame from wide format to long format, optionally
 |      leaving identifier variables set.
 |      
 |      This function is useful to massage a DataFrame into a format where one
 |      or more columns are identifier variables (`id_vars`), while all other
 |      columns, considered measured variables (`value_vars`), are "unpivoted" to
 |      the row axis, leaving just two non-identifier columns, 'variable' and
 |      'value'.
 |      
 |      .. versionadded:: 0.20.0
 |      
 |      Parameters
 |      ----------
 |      frame : DataFrame
 |      id_vars : tuple, list, or ndarray, optional
 |          Column(s) to use as identifier variables.
 |      value_vars : tuple, list, or ndarray, optional
 |          Column(s) to unpivot. If not specified, uses all columns that
 |          are not set as `id_vars`.
 |      var_name : scalar
 |          Name to use for the 'variable' column. If None it uses
 |          ``frame.columns.name`` or 'variable'.
 |      value_name : scalar, default 'value'
 |          Name to use for the 'value' column.
 |      col_level : int or string, optional
 |          If columns are a MultiIndex then use this level to melt.
 |      
 |      See also
 |      --------
 |      melt
 |      pivot_table
 |      DataFrame.pivot
 |      
 |      Examples
 |      --------
 |      >>> import pandas as pd
 |      >>> df = pd.DataFrame({'A': {0: 'a', 1: 'b', 2: 'c'},
 |      ...                    'B': {0: 1, 1: 3, 2: 5},
 |      ...                    'C': {0: 2, 1: 4, 2: 6}})
 |      >>> df
 |         A  B  C
 |      0  a  1  2
 |      1  b  3  4
 |      2  c  5  6
 |      
 |      >>> df.melt(id_vars=['A'], value_vars=['B'])
 |         A variable  value
 |      0  a        B      1
 |      1  b        B      3
 |      2  c        B      5
 |      
 |      >>> df.melt(id_vars=['A'], value_vars=['B', 'C'])
 |         A variable  value
 |      0  a        B      1
 |      1  b        B      3
 |      2  c        B      5
 |      3  a        C      2
 |      4  b        C      4
 |      5  c        C      6
 |      
 |      The names of 'variable' and 'value' columns can be customized:
 |      
 |      >>> df.melt(id_vars=['A'], value_vars=['B'],
 |      ...         var_name='myVarname', value_name='myValname')
 |         A myVarname  myValname
 |      0  a         B          1
 |      1  b         B          3
 |      2  c         B          5
 |      
 |      If you have multi-index columns:
 |      
 |      >>> df.columns = [list('ABC'), list('DEF')]
 |      >>> df
 |         A  B  C
 |         D  E  F
 |      0  a  1  2
 |      1  b  3  4
 |      2  c  5  6
 |      
 |      >>> df.melt(col_level=0, id_vars=['A'], value_vars=['B'])
 |         A variable  value
 |      0  a        B      1
 |      1  b        B      3
 |      2  c        B      5
 |      
 |      >>> df.melt(id_vars=[('A', 'D')], value_vars=[('B', 'E')])
 |        (A, D) variable_0 variable_1  value
 |      0      a          B          E      1
 |      1      b          B          E      3
 |      2      c          B          E      5
 |  
 |  memory_usage(self, index=True, deep=False)
 |      Return the memory usage of each column in bytes.
 |      
 |      The memory usage can optionally include the contribution of
 |      the index and elements of `object` dtype.
 |      
 |      This value is displayed in `DataFrame.info` by default. This can be
 |      suppressed by setting ``pandas.options.display.memory_usage`` to False.
 |      
 |      Parameters
 |      ----------
 |      index : bool, default True
 |          Specifies whether to include the memory usage of the DataFrame's
 |          index in returned Series. If ``index=True`` the memory usage of the
 |          index the first item in the output.
 |      deep : bool, default False
 |          If True, introspect the data deeply by interrogating
 |          `object` dtypes for system-level memory consumption, and include
 |          it in the returned values.
 |      
 |      Returns
 |      -------
 |      sizes : Series
 |          A Series whose index is the original column names and whose values
 |          is the memory usage of each column in bytes.
 |      
 |      See Also
 |      --------
 |      numpy.ndarray.nbytes : Total bytes consumed by the elements of an
 |          ndarray.
 |      Series.memory_usage : Bytes consumed by a Series.
 |      pandas.Categorical : Memory-efficient array for string values with
 |          many repeated values.
 |      DataFrame.info : Concise summary of a DataFrame.
 |      
 |      Examples
 |      --------
 |      >>> dtypes = ['int64', 'float64', 'complex128', 'object', 'bool']
 |      >>> data = dict([(t, np.ones(shape=5000).astype(t))
 |      ...              for t in dtypes])
 |      >>> df = pd.DataFrame(data)
 |      >>> df.head()
 |         int64  float64  complex128 object  bool
 |      0      1      1.0      (1+0j)      1  True
 |      1      1      1.0      (1+0j)      1  True
 |      2      1      1.0      (1+0j)      1  True
 |      3      1      1.0      (1+0j)      1  True
 |      4      1      1.0      (1+0j)      1  True
 |      
 |      >>> df.memory_usage()
 |      Index            80
 |      int64         40000
 |      float64       40000
 |      complex128    80000
 |      object        40000
 |      bool           5000
 |      dtype: int64
 |      
 |      >>> df.memory_usage(index=False)
 |      int64         40000
 |      float64       40000
 |      complex128    80000
 |      object        40000
 |      bool           5000
 |      dtype: int64
 |      
 |      The memory footprint of `object` dtype columns is ignored by default:
 |      
 |      >>> df.memory_usage(deep=True)
 |      Index             80
 |      int64          40000
 |      float64        40000
 |      complex128     80000
 |      object        160000
 |      bool            5000
 |      dtype: int64
 |      
 |      Use a Categorical for efficient storage of an object-dtype column with
 |      many repeated values.
 |      
 |      >>> df['object'].astype('category').memory_usage(deep=True)
 |      5168
 |  
 |  merge(self, right, how='inner', on=None, left_on=None, right_on=None, left_index=False, right_index=False, sort=False, suffixes=('_x', '_y'), copy=True, indicator=False, validate=None)
 |      Merge DataFrame objects by performing a database-style join operation by
 |      columns or indexes.
 |      
 |      If joining columns on columns, the DataFrame indexes *will be
 |      ignored*. Otherwise if joining indexes on indexes or indexes on a column or
 |      columns, the index will be passed on.
 |      
 |      Parameters
 |      ----------
 |      right : DataFrame
 |      how : {'left', 'right', 'outer', 'inner'}, default 'inner'
 |          * left: use only keys from left frame, similar to a SQL left outer join;
 |            preserve key order
 |          * right: use only keys from right frame, similar to a SQL right outer join;
 |            preserve key order
 |          * outer: use union of keys from both frames, similar to a SQL full outer
 |            join; sort keys lexicographically
 |          * inner: use intersection of keys from both frames, similar to a SQL inner
 |            join; preserve the order of the left keys
 |      on : label or list
 |          Column or index level names to join on. These must be found in both
 |          DataFrames. If `on` is None and not merging on indexes then this defaults
 |          to the intersection of the columns in both DataFrames.
 |      left_on : label or list, or array-like
 |          Column or index level names to join on in the left DataFrame. Can also
 |          be an array or list of arrays of the length of the left DataFrame.
 |          These arrays are treated as if they are columns.
 |      right_on : label or list, or array-like
 |          Column or index level names to join on in the right DataFrame. Can also
 |          be an array or list of arrays of the length of the right DataFrame.
 |          These arrays are treated as if they are columns.
 |      left_index : boolean, default False
 |          Use the index from the left DataFrame as the join key(s). If it is a
 |          MultiIndex, the number of keys in the other DataFrame (either the index
 |          or a number of columns) must match the number of levels
 |      right_index : boolean, default False
 |          Use the index from the right DataFrame as the join key. Same caveats as
 |          left_index
 |      sort : boolean, default False
 |          Sort the join keys lexicographically in the result DataFrame. If False,
 |          the order of the join keys depends on the join type (how keyword)
 |      suffixes : 2-length sequence (tuple, list, ...)
 |          Suffix to apply to overlapping column names in the left and right
 |          side, respectively
 |      copy : boolean, default True
 |          If False, do not copy data unnecessarily
 |      indicator : boolean or string, default False
 |          If True, adds a column to output DataFrame called "_merge" with
 |          information on the source of each row.
 |          If string, column with information on source of each row will be added to
 |          output DataFrame, and column will be named value of string.
 |          Information column is Categorical-type and takes on a value of "left_only"
 |          for observations whose merge key only appears in 'left' DataFrame,
 |          "right_only" for observations whose merge key only appears in 'right'
 |          DataFrame, and "both" if the observation's merge key is found in both.
 |      
 |      validate : string, default None
 |          If specified, checks if merge is of specified type.
 |      
 |          * "one_to_one" or "1:1": check if merge keys are unique in both
 |            left and right datasets.
 |          * "one_to_many" or "1:m": check if merge keys are unique in left
 |            dataset.
 |          * "many_to_one" or "m:1": check if merge keys are unique in right
 |            dataset.
 |          * "many_to_many" or "m:m": allowed, but does not result in checks.
 |      
 |          .. versionadded:: 0.21.0
 |      
 |      Notes
 |      -----
 |      Support for specifying index levels as the `on`, `left_on`, and
 |      `right_on` parameters was added in version 0.23.0
 |      
 |      Examples
 |      --------
 |      
 |      >>> A              >>> B
 |          lkey value         rkey value
 |      0   foo  1         0   foo  5
 |      1   bar  2         1   bar  6
 |      2   baz  3         2   qux  7
 |      3   foo  4         3   bar  8
 |      
 |      >>> A.merge(B, left_on='lkey', right_on='rkey', how='outer')
 |         lkey  value_x  rkey  value_y
 |      0  foo   1        foo   5
 |      1  foo   4        foo   5
 |      2  bar   2        bar   6
 |      3  bar   2        bar   8
 |      4  baz   3        NaN   NaN
 |      5  NaN   NaN      qux   7
 |      
 |      Returns
 |      -------
 |      merged : DataFrame
 |          The output type will the be same as 'left', if it is a subclass
 |          of DataFrame.
 |      
 |      See also
 |      --------
 |      merge_ordered
 |      merge_asof
 |      DataFrame.join
 |  
 |  min(self, axis=None, skipna=None, level=None, numeric_only=None, **kwargs)
 |      This method returns the minimum of the values in the object.
 |                  If you want the *index* of the minimum, use ``idxmin``. This is
 |                  the equivalent of the ``numpy.ndarray`` method ``argmin``.
 |      
 |      Parameters
 |      ----------
 |      axis : {index (0), columns (1)}
 |      skipna : boolean, default True
 |          Exclude NA/null values when computing the result.
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series
 |      numeric_only : boolean, default None
 |          Include only float, int, boolean columns. If None, will attempt to use
 |          everything, then use only numeric data. Not implemented for Series.
 |      
 |      Returns
 |      -------
 |      min : Series or DataFrame (if level specified)
 |  
 |  mod(self, other, axis='columns', level=None, fill_value=None)
 |      Modulo of dataframe and other, element-wise (binary operator `mod`).
 |      
 |      Equivalent to ``dataframe % other``, but with support to substitute a fill_value for
 |      missing data in one of the inputs.
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      None
 |      
 |      See also
 |      --------
 |      DataFrame.rmod
 |  
 |  mode(self, axis=0, numeric_only=False)
 |      Gets the mode(s) of each element along the axis selected. Adds a row
 |      for each mode per label, fills in gaps with nan.
 |      
 |      Note that there could be multiple values returned for the selected
 |      axis (when more than one item share the maximum frequency), which is
 |      the reason why a dataframe is returned. If you want to impute missing
 |      values with the mode in a dataframe ``df``, you can just do this:
 |      ``df.fillna(df.mode().iloc[0])``
 |      
 |      Parameters
 |      ----------
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          * 0 or 'index' : get mode of each column
 |          * 1 or 'columns' : get mode of each row
 |      numeric_only : boolean, default False
 |          if True, only apply to numeric columns
 |      
 |      Returns
 |      -------
 |      modes : DataFrame (sorted)
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'A': [1, 2, 1, 2, 1, 2, 3]})
 |      >>> df.mode()
 |         A
 |      0  1
 |      1  2
 |  
 |  mul(self, other, axis='columns', level=None, fill_value=None)
 |      Multiplication of dataframe and other, element-wise (binary operator `mul`).
 |      
 |      Equivalent to ``dataframe * other``, but with support to substitute a fill_value for
 |      missing data in one of the inputs.
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      None
 |      
 |      See also
 |      --------
 |      DataFrame.rmul
 |  
 |  multiply = mul(self, other, axis='columns', level=None, fill_value=None)
 |  
 |  ne(self, other, axis='columns', level=None)
 |      Wrapper for flexible comparison methods ne
 |  
 |  nlargest(self, n, columns, keep='first')
 |      Return the first `n` rows ordered by `columns` in descending order.
 |      
 |      Return the first `n` rows with the largest values in `columns`, in
 |      descending order. The columns that are not specified are returned as
 |      well, but not used for ordering.
 |      
 |      This method is equivalent to
 |      ``df.sort_values(columns, ascending=False).head(n)``, but more
 |      performant.
 |      
 |      Parameters
 |      ----------
 |      n : int
 |          Number of rows to return.
 |      columns : label or list of labels
 |          Column label(s) to order by.
 |      keep : {'first', 'last'}, default 'first'
 |          Where there are duplicate values:
 |      
 |          - `first` : prioritize the first occurrence(s)
 |          - `last` : prioritize the last occurrence(s)
 |      
 |      Returns
 |      -------
 |      DataFrame
 |          The first `n` rows ordered by the given columns in descending
 |          order.
 |      
 |      See Also
 |      --------
 |      DataFrame.nsmallest : Return the first `n` rows ordered by `columns` in
 |          ascending order.
 |      DataFrame.sort_values : Sort DataFrame by the values
 |      DataFrame.head : Return the first `n` rows without re-ordering.
 |      
 |      Notes
 |      -----
 |      This function cannot be used with all column types. For example, when
 |      specifying columns with `object` or `category` dtypes, ``TypeError`` is
 |      raised.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'a': [1, 10, 8, 10, -1],
 |      ...                    'b': list('abdce'),
 |      ...                    'c': [1.0, 2.0, np.nan, 3.0, 4.0]})
 |      >>> df
 |          a  b    c
 |      0   1  a  1.0
 |      1  10  b  2.0
 |      2   8  d  NaN
 |      3  10  c  3.0
 |      4  -1  e  4.0
 |      
 |      In the following example, we will use ``nlargest`` to select the three
 |      rows having the largest values in column "a".
 |      
 |      >>> df.nlargest(3, 'a')
 |          a  b    c
 |      1  10  b  2.0
 |      3  10  c  3.0
 |      2   8  d  NaN
 |      
 |      When using ``keep='last'``, ties are resolved in reverse order:
 |      
 |      >>> df.nlargest(3, 'a', keep='last')
 |          a  b    c
 |      3  10  c  3.0
 |      1  10  b  2.0
 |      2   8  d  NaN
 |      
 |      To order by the largest values in column "a" and then "c", we can
 |      specify multiple columns like in the next example.
 |      
 |      >>> df.nlargest(3, ['a', 'c'])
 |          a  b    c
 |      3  10  c  3.0
 |      1  10  b  2.0
 |      2   8  d  NaN
 |      
 |      Attempting to use ``nlargest`` on non-numeric dtypes will raise a
 |      ``TypeError``:
 |      
 |      >>> df.nlargest(3, 'b')
 |      Traceback (most recent call last):
 |      TypeError: Column 'b' has dtype object, cannot use method 'nlargest'
 |  
 |  notna(self)
 |      Detect existing (non-missing) values.
 |      
 |      Return a boolean same-sized object indicating if the values are not NA.
 |      Non-missing values get mapped to True. Characters such as empty
 |      strings ``''`` or :attr:`numpy.inf` are not considered NA values
 |      (unless you set ``pandas.options.mode.use_inf_as_na = True``).
 |      NA values, such as None or :attr:`numpy.NaN`, get mapped to False
 |      values.
 |      
 |      Returns
 |      -------
 |      DataFrame
 |          Mask of bool values for each element in DataFrame that
 |          indicates whether an element is not an NA value.
 |      
 |      See Also
 |      --------
 |      DataFrame.notnull : alias of notna
 |      DataFrame.isna : boolean inverse of notna
 |      DataFrame.dropna : omit axes labels with missing values
 |      notna : top-level notna
 |      
 |      Examples
 |      --------
 |      Show which entries in a DataFrame are not NA.
 |      
 |      >>> df = pd.DataFrame({'age': [5, 6, np.NaN],
 |      ...                    'born': [pd.NaT, pd.Timestamp('1939-05-27'),
 |      ...                             pd.Timestamp('1940-04-25')],
 |      ...                    'name': ['Alfred', 'Batman', ''],
 |      ...                    'toy': [None, 'Batmobile', 'Joker']})
 |      >>> df
 |         age       born    name        toy
 |      0  5.0        NaT  Alfred       None
 |      1  6.0 1939-05-27  Batman  Batmobile
 |      2  NaN 1940-04-25              Joker
 |      
 |      >>> df.notna()
 |           age   born  name    toy
 |      0   True  False  True  False
 |      1   True   True  True   True
 |      2  False   True  True   True
 |      
 |      Show which entries in a Series are not NA.
 |      
 |      >>> ser = pd.Series([5, 6, np.NaN])
 |      >>> ser
 |      0    5.0
 |      1    6.0
 |      2    NaN
 |      dtype: float64
 |      
 |      >>> ser.notna()
 |      0     True
 |      1     True
 |      2    False
 |      dtype: bool
 |  
 |  notnull(self)
 |      Detect existing (non-missing) values.
 |      
 |      Return a boolean same-sized object indicating if the values are not NA.
 |      Non-missing values get mapped to True. Characters such as empty
 |      strings ``''`` or :attr:`numpy.inf` are not considered NA values
 |      (unless you set ``pandas.options.mode.use_inf_as_na = True``).
 |      NA values, such as None or :attr:`numpy.NaN`, get mapped to False
 |      values.
 |      
 |      Returns
 |      -------
 |      DataFrame
 |          Mask of bool values for each element in DataFrame that
 |          indicates whether an element is not an NA value.
 |      
 |      See Also
 |      --------
 |      DataFrame.notnull : alias of notna
 |      DataFrame.isna : boolean inverse of notna
 |      DataFrame.dropna : omit axes labels with missing values
 |      notna : top-level notna
 |      
 |      Examples
 |      --------
 |      Show which entries in a DataFrame are not NA.
 |      
 |      >>> df = pd.DataFrame({'age': [5, 6, np.NaN],
 |      ...                    'born': [pd.NaT, pd.Timestamp('1939-05-27'),
 |      ...                             pd.Timestamp('1940-04-25')],
 |      ...                    'name': ['Alfred', 'Batman', ''],
 |      ...                    'toy': [None, 'Batmobile', 'Joker']})
 |      >>> df
 |         age       born    name        toy
 |      0  5.0        NaT  Alfred       None
 |      1  6.0 1939-05-27  Batman  Batmobile
 |      2  NaN 1940-04-25              Joker
 |      
 |      >>> df.notna()
 |           age   born  name    toy
 |      0   True  False  True  False
 |      1   True   True  True   True
 |      2  False   True  True   True
 |      
 |      Show which entries in a Series are not NA.
 |      
 |      >>> ser = pd.Series([5, 6, np.NaN])
 |      >>> ser
 |      0    5.0
 |      1    6.0
 |      2    NaN
 |      dtype: float64
 |      
 |      >>> ser.notna()
 |      0     True
 |      1     True
 |      2    False
 |      dtype: bool
 |  
 |  nsmallest(self, n, columns, keep='first')
 |      Get the rows of a DataFrame sorted by the `n` smallest
 |      values of `columns`.
 |      
 |      Parameters
 |      ----------
 |      n : int
 |          Number of items to retrieve
 |      columns : list or str
 |          Column name or names to order by
 |      keep : {'first', 'last'}, default 'first'
 |          Where there are duplicate values:
 |          - ``first`` : take the first occurrence.
 |          - ``last`` : take the last occurrence.
 |      
 |      Returns
 |      -------
 |      DataFrame
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'a': [1, 10, 8, 11, -1],
 |      ...                    'b': list('abdce'),
 |      ...                    'c': [1.0, 2.0, np.nan, 3.0, 4.0]})
 |      >>> df.nsmallest(3, 'a')
 |         a  b   c
 |      4 -1  e   4
 |      0  1  a   1
 |      2  8  d NaN
 |  
 |  nunique(self, axis=0, dropna=True)
 |      Return Series with number of distinct observations over requested
 |      axis.
 |      
 |      .. versionadded:: 0.20.0
 |      
 |      Parameters
 |      ----------
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |      dropna : boolean, default True
 |          Don't include NaN in the counts.
 |      
 |      Returns
 |      -------
 |      nunique : Series
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'A': [1, 2, 3], 'B': [1, 1, 1]})
 |      >>> df.nunique()
 |      A    3
 |      B    1
 |      
 |      >>> df.nunique(axis=1)
 |      0    1
 |      1    2
 |      2    2
 |  
 |  pivot(self, index=None, columns=None, values=None)
 |      Return reshaped DataFrame organized by given index / column values.
 |      
 |      Reshape data (produce a "pivot" table) based on column values. Uses
 |      unique values from specified `index` / `columns` to form axes of the
 |      resulting DataFrame. This function does not support data
 |      aggregation, multiple values will result in a MultiIndex in the
 |      columns. See the :ref:`User Guide <reshaping>` for more on reshaping.
 |      
 |      Parameters
 |      ----------
 |      index : string or object, optional
 |          Column to use to make new frame's index. If None, uses
 |          existing index.
 |      columns : string or object
 |          Column to use to make new frame's columns.
 |      values : string, object or a list of the previous, optional
 |          Column(s) to use for populating new frame's values. If not
 |          specified, all remaining columns will be used and the result will
 |          have hierarchically indexed columns.
 |      
 |          .. versionchanged :: 0.23.0
 |             Also accept list of column names.
 |      
 |      Returns
 |      -------
 |      DataFrame
 |          Returns reshaped DataFrame.
 |      
 |      Raises
 |      ------
 |      ValueError:
 |          When there are any `index`, `columns` combinations with multiple
 |          values. `DataFrame.pivot_table` when you need to aggregate.
 |      
 |      See Also
 |      --------
 |      DataFrame.pivot_table : generalization of pivot that can handle
 |          duplicate values for one index/column pair.
 |      DataFrame.unstack : pivot based on the index values instead of a
 |          column.
 |      
 |      Notes
 |      -----
 |      For finer-tuned control, see hierarchical indexing documentation along
 |      with the related stack/unstack methods.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'foo': ['one', 'one', 'one', 'two', 'two',
 |      ...                            'two'],
 |      ...                    'bar': ['A', 'B', 'C', 'A', 'B', 'C'],
 |      ...                    'baz': [1, 2, 3, 4, 5, 6],
 |      ...                    'zoo': ['x', 'y', 'z', 'q', 'w', 't']})
 |      >>> df
 |          foo   bar  baz  zoo
 |      0   one   A    1    x
 |      1   one   B    2    y
 |      2   one   C    3    z
 |      3   two   A    4    q
 |      4   two   B    5    w
 |      5   two   C    6    t
 |      
 |      >>> df.pivot(index='foo', columns='bar', values='baz')
 |      bar  A   B   C
 |      foo
 |      one  1   2   3
 |      two  4   5   6
 |      
 |      >>> df.pivot(index='foo', columns='bar')['baz']
 |      bar  A   B   C
 |      foo
 |      one  1   2   3
 |      two  4   5   6
 |      
 |      >>> df.pivot(index='foo', columns='bar', values=['baz', 'zoo'])
 |            baz       zoo
 |      bar   A  B  C   A  B  C
 |      foo
 |      one   1  2  3   x  y  z
 |      two   4  5  6   q  w  t
 |      
 |      A ValueError is raised if there are any duplicates.
 |      
 |      >>> df = pd.DataFrame({"foo": ['one', 'one', 'two', 'two'],
 |      ...                    "bar": ['A', 'A', 'B', 'C'],
 |      ...                    "baz": [1, 2, 3, 4]})
 |      >>> df
 |         foo bar  baz
 |      0  one   A    1
 |      1  one   A    2
 |      2  two   B    3
 |      3  two   C    4
 |      
 |      Notice that the first two rows are the same for our `index`
 |      and `columns` arguments.
 |      
 |      >>> df.pivot(index='foo', columns='bar', values='baz')
 |      Traceback (most recent call last):
 |         ...
 |      ValueError: Index contains duplicate entries, cannot reshape
 |  
 |  pivot_table(self, values=None, index=None, columns=None, aggfunc='mean', fill_value=None, margins=False, dropna=True, margins_name='All')
 |      Create a spreadsheet-style pivot table as a DataFrame. The levels in
 |      the pivot table will be stored in MultiIndex objects (hierarchical
 |      indexes) on the index and columns of the result DataFrame
 |      
 |      Parameters
 |      ----------
 |      values : column to aggregate, optional
 |      index : column, Grouper, array, or list of the previous
 |          If an array is passed, it must be the same length as the data. The
 |          list can contain any of the other types (except list).
 |          Keys to group by on the pivot table index.  If an array is passed,
 |          it is being used as the same manner as column values.
 |      columns : column, Grouper, array, or list of the previous
 |          If an array is passed, it must be the same length as the data. The
 |          list can contain any of the other types (except list).
 |          Keys to group by on the pivot table column.  If an array is passed,
 |          it is being used as the same manner as column values.
 |      aggfunc : function, list of functions, dict, default numpy.mean
 |          If list of functions passed, the resulting pivot table will have
 |          hierarchical columns whose top level are the function names
 |          (inferred from the function objects themselves)
 |          If dict is passed, the key is column to aggregate and value
 |          is function or list of functions
 |      fill_value : scalar, default None
 |          Value to replace missing values with
 |      margins : boolean, default False
 |          Add all row / columns (e.g. for subtotal / grand totals)
 |      dropna : boolean, default True
 |          Do not include columns whose entries are all NaN
 |      margins_name : string, default 'All'
 |          Name of the row / column that will contain the totals
 |          when margins is True.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({"A": ["foo", "foo", "foo", "foo", "foo",
 |      ...                          "bar", "bar", "bar", "bar"],
 |      ...                    "B": ["one", "one", "one", "two", "two",
 |      ...                          "one", "one", "two", "two"],
 |      ...                    "C": ["small", "large", "large", "small",
 |      ...                          "small", "large", "small", "small",
 |      ...                          "large"],
 |      ...                    "D": [1, 2, 2, 3, 3, 4, 5, 6, 7]})
 |      >>> df
 |           A    B      C  D
 |      0  foo  one  small  1
 |      1  foo  one  large  2
 |      2  foo  one  large  2
 |      3  foo  two  small  3
 |      4  foo  two  small  3
 |      5  bar  one  large  4
 |      6  bar  one  small  5
 |      7  bar  two  small  6
 |      8  bar  two  large  7
 |      
 |      >>> table = pivot_table(df, values='D', index=['A', 'B'],
 |      ...                     columns=['C'], aggfunc=np.sum)
 |      >>> table
 |      C        large  small
 |      A   B
 |      bar one    4.0    5.0
 |          two    7.0    6.0
 |      foo one    4.0    1.0
 |          two    NaN    6.0
 |      
 |      >>> table = pivot_table(df, values='D', index=['A', 'B'],
 |      ...                     columns=['C'], aggfunc=np.sum)
 |      >>> table
 |      C        large  small
 |      A   B
 |      bar one    4.0    5.0
 |          two    7.0    6.0
 |      foo one    4.0    1.0
 |          two    NaN    6.0
 |      
 |      >>> table = pivot_table(df, values=['D', 'E'], index=['A', 'C'],
 |      ...                     aggfunc={'D': np.mean,
 |      ...                              'E': [min, max, np.mean]})
 |      >>> table
 |                        D   E
 |                     mean max median min
 |      A   C
 |      bar large  5.500000  16   14.5  13
 |          small  5.500000  15   14.5  14
 |      foo large  2.000000  10    9.5   9
 |          small  2.333333  12   11.0   8
 |      
 |      Returns
 |      -------
 |      table : DataFrame
 |      
 |      See also
 |      --------
 |      DataFrame.pivot : pivot without aggregation that can handle
 |          non-numeric data
 |  
 |  pow(self, other, axis='columns', level=None, fill_value=None)
 |      Exponential power of dataframe and other, element-wise (binary operator `pow`).
 |      
 |      Equivalent to ``dataframe ** other``, but with support to substitute a fill_value for
 |      missing data in one of the inputs.
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      None
 |      
 |      See also
 |      --------
 |      DataFrame.rpow
 |  
 |  prod(self, axis=None, skipna=None, level=None, numeric_only=None, min_count=0, **kwargs)
 |      Return the product of the values for the requested axis
 |      
 |      Parameters
 |      ----------
 |      axis : {index (0), columns (1)}
 |      skipna : boolean, default True
 |          Exclude NA/null values when computing the result.
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series
 |      numeric_only : boolean, default None
 |          Include only float, int, boolean columns. If None, will attempt to use
 |          everything, then use only numeric data. Not implemented for Series.
 |      min_count : int, default 0
 |          The required number of valid values to perform the operation. If fewer than
 |          ``min_count`` non-NA values are present the result will be NA.
 |      
 |          .. versionadded :: 0.22.0
 |      
 |             Added with the default being 0. This means the sum of an all-NA
 |             or empty Series is 0, and the product of an all-NA or empty
 |             Series is 1.
 |      
 |      Returns
 |      -------
 |      prod : Series or DataFrame (if level specified)
 |      
 |      Examples
 |      --------
 |      By default, the product of an empty or all-NA Series is ``1``
 |      
 |      >>> pd.Series([]).prod()
 |      1.0
 |      
 |      This can be controlled with the ``min_count`` parameter
 |      
 |      >>> pd.Series([]).prod(min_count=1)
 |      nan
 |      
 |      Thanks to the ``skipna`` parameter, ``min_count`` handles all-NA and
 |      empty series identically.
 |      
 |      >>> pd.Series([np.nan]).prod()
 |      1.0
 |      
 |      >>> pd.Series([np.nan]).prod(min_count=1)
 |      nan
 |  
 |  product = prod(self, axis=None, skipna=None, level=None, numeric_only=None, min_count=0, **kwargs)
 |  
 |  quantile(self, q=0.5, axis=0, numeric_only=True, interpolation='linear')
 |      Return values at the given quantile over requested axis, a la
 |      numpy.percentile.
 |      
 |      Parameters
 |      ----------
 |      q : float or array-like, default 0.5 (50% quantile)
 |          0 <= q <= 1, the quantile(s) to compute
 |      axis : {0, 1, 'index', 'columns'} (default 0)
 |          0 or 'index' for row-wise, 1 or 'columns' for column-wise
 |      numeric_only : boolean, default True
 |          If False, the quantile of datetime and timedelta data will be
 |          computed as well
 |      interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
 |          .. versionadded:: 0.18.0
 |      
 |          This optional parameter specifies the interpolation method to use,
 |          when the desired quantile lies between two data points `i` and `j`:
 |      
 |          * linear: `i + (j - i) * fraction`, where `fraction` is the
 |            fractional part of the index surrounded by `i` and `j`.
 |          * lower: `i`.
 |          * higher: `j`.
 |          * nearest: `i` or `j` whichever is nearest.
 |          * midpoint: (`i` + `j`) / 2.
 |      
 |      Returns
 |      -------
 |      quantiles : Series or DataFrame
 |      
 |          - If ``q`` is an array, a DataFrame will be returned where the
 |            index is ``q``, the columns are the columns of self, and the
 |            values are the quantiles.
 |          - If ``q`` is a float, a Series will be returned where the
 |            index is the columns of self and the values are the quantiles.
 |      
 |      Examples
 |      --------
 |      
 |      >>> df = pd.DataFrame(np.array([[1, 1], [2, 10], [3, 100], [4, 100]]),
 |                            columns=['a', 'b'])
 |      >>> df.quantile(.1)
 |      a    1.3
 |      b    3.7
 |      dtype: float64
 |      >>> df.quantile([.1, .5])
 |             a     b
 |      0.1  1.3   3.7
 |      0.5  2.5  55.0
 |      
 |      Specifying `numeric_only=False` will also compute the quantile of
 |      datetime and timedelta data.
 |      
 |      >>> df = pd.DataFrame({'A': [1, 2],
 |                             'B': [pd.Timestamp('2010'),
 |                                   pd.Timestamp('2011')],
 |                             'C': [pd.Timedelta('1 days'),
 |                                   pd.Timedelta('2 days')]})
 |      >>> df.quantile(0.5, numeric_only=False)
 |      A                    1.5
 |      B    2010-07-02 12:00:00
 |      C        1 days 12:00:00
 |      Name: 0.5, dtype: object
 |      
 |      See Also
 |      --------
 |      pandas.core.window.Rolling.quantile
 |  
 |  query(self, expr, inplace=False, **kwargs)
 |      Query the columns of a frame with a boolean expression.
 |      
 |      Parameters
 |      ----------
 |      expr : string
 |          The query string to evaluate.  You can refer to variables
 |          in the environment by prefixing them with an '@' character like
 |          ``@a + b``.
 |      inplace : bool
 |          Whether the query should modify the data in place or return
 |          a modified copy
 |      
 |          .. versionadded:: 0.18.0
 |      
 |      kwargs : dict
 |          See the documentation for :func:`pandas.eval` for complete details
 |          on the keyword arguments accepted by :meth:`DataFrame.query`.
 |      
 |      Returns
 |      -------
 |      q : DataFrame
 |      
 |      Notes
 |      -----
 |      The result of the evaluation of this expression is first passed to
 |      :attr:`DataFrame.loc` and if that fails because of a
 |      multidimensional key (e.g., a DataFrame) then the result will be passed
 |      to :meth:`DataFrame.__getitem__`.
 |      
 |      This method uses the top-level :func:`pandas.eval` function to
 |      evaluate the passed query.
 |      
 |      The :meth:`~pandas.DataFrame.query` method uses a slightly
 |      modified Python syntax by default. For example, the ``&`` and ``|``
 |      (bitwise) operators have the precedence of their boolean cousins,
 |      :keyword:`and` and :keyword:`or`. This *is* syntactically valid Python,
 |      however the semantics are different.
 |      
 |      You can change the semantics of the expression by passing the keyword
 |      argument ``parser='python'``. This enforces the same semantics as
 |      evaluation in Python space. Likewise, you can pass ``engine='python'``
 |      to evaluate an expression using Python itself as a backend. This is not
 |      recommended as it is inefficient compared to using ``numexpr`` as the
 |      engine.
 |      
 |      The :attr:`DataFrame.index` and
 |      :attr:`DataFrame.columns` attributes of the
 |      :class:`~pandas.DataFrame` instance are placed in the query namespace
 |      by default, which allows you to treat both the index and columns of the
 |      frame as a column in the frame.
 |      The identifier ``index`` is used for the frame index; you can also
 |      use the name of the index to identify it in a query. Please note that
 |      Python keywords may not be used as identifiers.
 |      
 |      For further details and examples see the ``query`` documentation in
 |      :ref:`indexing <indexing.query>`.
 |      
 |      See Also
 |      --------
 |      pandas.eval
 |      DataFrame.eval
 |      
 |      Examples
 |      --------
 |      >>> from numpy.random import randn
 |      >>> from pandas import DataFrame
 |      >>> df = pd.DataFrame(randn(10, 2), columns=list('ab'))
 |      >>> df.query('a > b')
 |      >>> df[df.a > df.b]  # same result as the previous expression
 |  
 |  radd(self, other, axis='columns', level=None, fill_value=None)
 |      Addition of dataframe and other, element-wise (binary operator `radd`).
 |      
 |      Equivalent to ``other + dataframe``, but with support to substitute a fill_value for
 |      missing data in one of the inputs.
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      
 |      >>> a = pd.DataFrame([1, 1, 1, np.nan], index=['a', 'b', 'c', 'd'],
 |      ...                  columns=['one'])
 |      >>> a
 |         one
 |      a  1.0
 |      b  1.0
 |      c  1.0
 |      d  NaN
 |      >>> b = pd.DataFrame(dict(one=[1, np.nan, 1, np.nan],
 |      ...                       two=[np.nan, 2, np.nan, 2]),
 |      ...                  index=['a', 'b', 'd', 'e'])
 |      >>> b
 |         one  two
 |      a  1.0  NaN
 |      b  NaN  2.0
 |      d  1.0  NaN
 |      e  NaN  2.0
 |      >>> a.add(b, fill_value=0)
 |         one  two
 |      a  2.0  NaN
 |      b  1.0  2.0
 |      c  1.0  NaN
 |      d  1.0  NaN
 |      e  NaN  2.0
 |      
 |      
 |      See also
 |      --------
 |      DataFrame.add
 |  
 |  rdiv = rtruediv(self, other, axis='columns', level=None, fill_value=None)
 |  
 |  reindex(self, labels=None, index=None, columns=None, axis=None, method=None, copy=True, level=None, fill_value=nan, limit=None, tolerance=None)
 |      Conform DataFrame to new index with optional filling logic, placing
 |      NA/NaN in locations having no value in the previous index. A new object
 |      is produced unless the new index is equivalent to the current one and
 |      copy=False
 |      
 |      Parameters
 |      ----------
 |      labels : array-like, optional
 |          New labels / index to conform the axis specified by 'axis' to.
 |      index, columns : array-like, optional (should be specified using keywords)
 |          New labels / index to conform to. Preferably an Index object to
 |          avoid duplicating data
 |      axis : int or str, optional
 |          Axis to target. Can be either the axis name ('index', 'columns')
 |          or number (0, 1).
 |      method : {None, 'backfill'/'bfill', 'pad'/'ffill', 'nearest'}, optional
 |          method to use for filling holes in reindexed DataFrame.
 |          Please note: this is only applicable to DataFrames/Series with a
 |          monotonically increasing/decreasing index.
 |      
 |          * default: don't fill gaps
 |          * pad / ffill: propagate last valid observation forward to next
 |            valid
 |          * backfill / bfill: use next valid observation to fill gap
 |          * nearest: use nearest valid observations to fill gap
 |      
 |      copy : boolean, default True
 |          Return a new object, even if the passed indexes are the same
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : scalar, default np.NaN
 |          Value to use for missing values. Defaults to NaN, but can be any
 |          "compatible" value
 |      limit : int, default None
 |          Maximum number of consecutive elements to forward or backward fill
 |      tolerance : optional
 |          Maximum distance between original and new labels for inexact
 |          matches. The values of the index at the matching locations most
 |          satisfy the equation ``abs(index[indexer] - target) <= tolerance``.
 |      
 |          Tolerance may be a scalar value, which applies the same tolerance
 |          to all values, or list-like, which applies variable tolerance per
 |          element. List-like includes list, tuple, array, Series, and must be
 |          the same size as the index and its dtype must exactly match the
 |          index's type.
 |      
 |          .. versionadded:: 0.21.0 (list-like tolerance)
 |      
 |      Examples
 |      --------
 |      
 |      ``DataFrame.reindex`` supports two calling conventions
 |      
 |      * ``(index=index_labels, columns=column_labels, ...)``
 |      * ``(labels, axis={'index', 'columns'}, ...)``
 |      
 |      We *highly* recommend using keyword arguments to clarify your
 |      intent.
 |      
 |      Create a dataframe with some fictional data.
 |      
 |      >>> index = ['Firefox', 'Chrome', 'Safari', 'IE10', 'Konqueror']
 |      >>> df = pd.DataFrame({
 |      ...      'http_status': [200,200,404,404,301],
 |      ...      'response_time': [0.04, 0.02, 0.07, 0.08, 1.0]},
 |      ...       index=index)
 |      >>> df
 |                 http_status  response_time
 |      Firefox            200           0.04
 |      Chrome             200           0.02
 |      Safari             404           0.07
 |      IE10               404           0.08
 |      Konqueror          301           1.00
 |      
 |      Create a new index and reindex the dataframe. By default
 |      values in the new index that do not have corresponding
 |      records in the dataframe are assigned ``NaN``.
 |      
 |      >>> new_index= ['Safari', 'Iceweasel', 'Comodo Dragon', 'IE10',
 |      ...             'Chrome']
 |      >>> df.reindex(new_index)
 |                     http_status  response_time
 |      Safari               404.0           0.07
 |      Iceweasel              NaN            NaN
 |      Comodo Dragon          NaN            NaN
 |      IE10                 404.0           0.08
 |      Chrome               200.0           0.02
 |      
 |      We can fill in the missing values by passing a value to
 |      the keyword ``fill_value``. Because the index is not monotonically
 |      increasing or decreasing, we cannot use arguments to the keyword
 |      ``method`` to fill the ``NaN`` values.
 |      
 |      >>> df.reindex(new_index, fill_value=0)
 |                     http_status  response_time
 |      Safari                 404           0.07
 |      Iceweasel                0           0.00
 |      Comodo Dragon            0           0.00
 |      IE10                   404           0.08
 |      Chrome                 200           0.02
 |      
 |      >>> df.reindex(new_index, fill_value='missing')
 |                    http_status response_time
 |      Safari                404          0.07
 |      Iceweasel         missing       missing
 |      Comodo Dragon     missing       missing
 |      IE10                  404          0.08
 |      Chrome                200          0.02
 |      
 |      We can also reindex the columns.
 |      
 |      >>> df.reindex(columns=['http_status', 'user_agent'])
 |                 http_status  user_agent
 |      Firefox            200         NaN
 |      Chrome             200         NaN
 |      Safari             404         NaN
 |      IE10               404         NaN
 |      Konqueror          301         NaN
 |      
 |      Or we can use "axis-style" keyword arguments
 |      
 |      >>> df.reindex(['http_status', 'user_agent'], axis="columns")
 |                 http_status  user_agent
 |      Firefox            200         NaN
 |      Chrome             200         NaN
 |      Safari             404         NaN
 |      IE10               404         NaN
 |      Konqueror          301         NaN
 |      
 |      To further illustrate the filling functionality in
 |      ``reindex``, we will create a dataframe with a
 |      monotonically increasing index (for example, a sequence
 |      of dates).
 |      
 |      >>> date_index = pd.date_range('1/1/2010', periods=6, freq='D')
 |      >>> df2 = pd.DataFrame({"prices": [100, 101, np.nan, 100, 89, 88]},
 |      ...                    index=date_index)
 |      >>> df2
 |                  prices
 |      2010-01-01     100
 |      2010-01-02     101
 |      2010-01-03     NaN
 |      2010-01-04     100
 |      2010-01-05      89
 |      2010-01-06      88
 |      
 |      Suppose we decide to expand the dataframe to cover a wider
 |      date range.
 |      
 |      >>> date_index2 = pd.date_range('12/29/2009', periods=10, freq='D')
 |      >>> df2.reindex(date_index2)
 |                  prices
 |      2009-12-29     NaN
 |      2009-12-30     NaN
 |      2009-12-31     NaN
 |      2010-01-01     100
 |      2010-01-02     101
 |      2010-01-03     NaN
 |      2010-01-04     100
 |      2010-01-05      89
 |      2010-01-06      88
 |      2010-01-07     NaN
 |      
 |      The index entries that did not have a value in the original data frame
 |      (for example, '2009-12-29') are by default filled with ``NaN``.
 |      If desired, we can fill in the missing values using one of several
 |      options.
 |      
 |      For example, to backpropagate the last valid value to fill the ``NaN``
 |      values, pass ``bfill`` as an argument to the ``method`` keyword.
 |      
 |      >>> df2.reindex(date_index2, method='bfill')
 |                  prices
 |      2009-12-29     100
 |      2009-12-30     100
 |      2009-12-31     100
 |      2010-01-01     100
 |      2010-01-02     101
 |      2010-01-03     NaN
 |      2010-01-04     100
 |      2010-01-05      89
 |      2010-01-06      88
 |      2010-01-07     NaN
 |      
 |      Please note that the ``NaN`` value present in the original dataframe
 |      (at index value 2010-01-03) will not be filled by any of the
 |      value propagation schemes. This is because filling while reindexing
 |      does not look at dataframe values, but only compares the original and
 |      desired indexes. If you do want to fill in the ``NaN`` values present
 |      in the original dataframe, use the ``fillna()`` method.
 |      
 |      See the :ref:`user guide <basics.reindexing>` for more.
 |      
 |      Returns
 |      -------
 |      reindexed : DataFrame
 |  
 |  reindex_axis(self, labels, axis=0, method=None, level=None, copy=True, limit=None, fill_value=nan)
 |      Conform input object to new index with optional
 |      filling logic, placing NA/NaN in locations having no value in the
 |      previous index. A new object is produced unless the new index is
 |      equivalent to the current one and copy=False
 |      
 |      Parameters
 |      ----------
 |      labels : array-like
 |          New labels / index to conform to. Preferably an Index object to
 |          avoid duplicating data
 |      axis : {0 or 'index', 1 or 'columns'}
 |      method : {None, 'backfill'/'bfill', 'pad'/'ffill', 'nearest'}, optional
 |          Method to use for filling holes in reindexed DataFrame:
 |      
 |          * default: don't fill gaps
 |          * pad / ffill: propagate last valid observation forward to next
 |            valid
 |          * backfill / bfill: use next valid observation to fill gap
 |          * nearest: use nearest valid observations to fill gap
 |      
 |      copy : boolean, default True
 |          Return a new object, even if the passed indexes are the same
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      limit : int, default None
 |          Maximum number of consecutive elements to forward or backward fill
 |      tolerance : optional
 |          Maximum distance between original and new labels for inexact
 |          matches. The values of the index at the matching locations most
 |          satisfy the equation ``abs(index[indexer] - target) <= tolerance``.
 |      
 |          Tolerance may be a scalar value, which applies the same tolerance
 |          to all values, or list-like, which applies variable tolerance per
 |          element. List-like includes list, tuple, array, Series, and must be
 |          the same size as the index and its dtype must exactly match the
 |          index's type.
 |      
 |          .. versionadded:: 0.21.0 (list-like tolerance)
 |      
 |      Examples
 |      --------
 |      >>> df.reindex_axis(['A', 'B', 'C'], axis=1)
 |      
 |      See Also
 |      --------
 |      reindex, reindex_like
 |      
 |      Returns
 |      -------
 |      reindexed : DataFrame
 |  
 |  rename(self, mapper=None, index=None, columns=None, axis=None, copy=True, inplace=False, level=None)
 |      Alter axes labels.
 |      
 |      Function / dict values must be unique (1-to-1). Labels not contained in
 |      a dict / Series will be left as-is. Extra labels listed don't throw an
 |      error.
 |      
 |      See the :ref:`user guide <basics.rename>` for more.
 |      
 |      Parameters
 |      ----------
 |      mapper, index, columns : dict-like or function, optional
 |          dict-like or functions transformations to apply to
 |          that axis' values. Use either ``mapper`` and ``axis`` to
 |          specify the axis to target with ``mapper``, or ``index`` and
 |          ``columns``.
 |      axis : int or str, optional
 |          Axis to target with ``mapper``. Can be either the axis name
 |          ('index', 'columns') or number (0, 1). The default is 'index'.
 |      copy : boolean, default True
 |          Also copy underlying data
 |      inplace : boolean, default False
 |          Whether to return a new DataFrame. If True then value of copy is
 |          ignored.
 |      level : int or level name, default None
 |          In case of a MultiIndex, only rename labels in the specified
 |          level.
 |      
 |      Returns
 |      -------
 |      renamed : DataFrame
 |      
 |      See Also
 |      --------
 |      pandas.DataFrame.rename_axis
 |      
 |      Examples
 |      --------
 |      
 |      ``DataFrame.rename`` supports two calling conventions
 |      
 |      * ``(index=index_mapper, columns=columns_mapper, ...)``
 |      * ``(mapper, axis={'index', 'columns'}, ...)``
 |      
 |      We *highly* recommend using keyword arguments to clarify your
 |      intent.
 |      
 |      >>> df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
 |      >>> df.rename(index=str, columns={"A": "a", "B": "c"})
 |         a  c
 |      0  1  4
 |      1  2  5
 |      2  3  6
 |      
 |      >>> df.rename(index=str, columns={"A": "a", "C": "c"})
 |         a  B
 |      0  1  4
 |      1  2  5
 |      2  3  6
 |      
 |      Using axis-style parameters
 |      
 |      >>> df.rename(str.lower, axis='columns')
 |         a  b
 |      0  1  4
 |      1  2  5
 |      2  3  6
 |      
 |      >>> df.rename({1: 2, 2: 4}, axis='index')
 |         A  B
 |      0  1  4
 |      2  2  5
 |      4  3  6
 |  
 |  reorder_levels(self, order, axis=0)
 |      Rearrange index levels using input order.
 |      May not drop or duplicate levels
 |      
 |      Parameters
 |      ----------
 |      order : list of int or list of str
 |          List representing new level order. Reference level by number
 |          (position) or by key (label).
 |      axis : int
 |          Where to reorder levels.
 |      
 |      Returns
 |      -------
 |      type of caller (new object)
 |  
 |  replace(self, to_replace=None, value=None, inplace=False, limit=None, regex=False, method='pad')
 |      Replace values given in `to_replace` with `value`.
 |      
 |      Values of the DataFrame are replaced with other values dynamically.
 |      This differs from updating with ``.loc`` or ``.iloc``, which require
 |      you to specify a location to update with some value.
 |      
 |      Parameters
 |      ----------
 |      to_replace : str, regex, list, dict, Series, int, float, or None
 |          How to find the values that will be replaced.
 |      
 |          * numeric, str or regex:
 |      
 |              - numeric: numeric values equal to `to_replace` will be
 |                replaced with `value`
 |              - str: string exactly matching `to_replace` will be replaced
 |                with `value`
 |              - regex: regexs matching `to_replace` will be replaced with
 |                `value`
 |      
 |          * list of str, regex, or numeric:
 |      
 |              - First, if `to_replace` and `value` are both lists, they
 |                **must** be the same length.
 |              - Second, if ``regex=True`` then all of the strings in **both**
 |                lists will be interpreted as regexs otherwise they will match
 |                directly. This doesn't matter much for `value` since there
 |                are only a few possible substitution regexes you can use.
 |              - str, regex and numeric rules apply as above.
 |      
 |          * dict:
 |      
 |              - Dicts can be used to specify different replacement values
 |                for different existing values. For example,
 |                ``{'a': 'b', 'y': 'z'}`` replaces the value 'a' with 'b' and
 |                'y' with 'z'. To use a dict in this way the `value`
 |                parameter should be `None`.
 |              - For a DataFrame a dict can specify that different values
 |                should be replaced in different columns. For example,
 |                ``{'a': 1, 'b': 'z'}`` looks for the value 1 in column 'a'
 |                and the value 'z' in column 'b' and replaces these values
 |                with whatever is specified in `value`. The `value` parameter
 |                should not be ``None`` in this case. You can treat this as a
 |                special case of passing two lists except that you are
 |                specifying the column to search in.
 |              - For a DataFrame nested dictionaries, e.g.,
 |                ``{'a': {'b': np.nan}}``, are read as follows: look in column
 |                'a' for the value 'b' and replace it with NaN. The `value`
 |                parameter should be ``None`` to use a nested dict in this
 |                way. You can nest regular expressions as well. Note that
 |                column names (the top-level dictionary keys in a nested
 |                dictionary) **cannot** be regular expressions.
 |      
 |          * None:
 |      
 |              - This means that the `regex` argument must be a string,
 |                compiled regular expression, or list, dict, ndarray or
 |                Series of such elements. If `value` is also ``None`` then
 |                this **must** be a nested dictionary or Series.
 |      
 |          See the examples section for examples of each of these.
 |      value : scalar, dict, list, str, regex, default None
 |          Value to replace any values matching `to_replace` with.
 |          For a DataFrame a dict of values can be used to specify which
 |          value to use for each column (columns not in the dict will not be
 |          filled). Regular expressions, strings and lists or dicts of such
 |          objects are also allowed.
 |      inplace : boolean, default False
 |          If True, in place. Note: this will modify any
 |          other views on this object (e.g. a column from a DataFrame).
 |          Returns the caller if this is True.
 |      limit : int, default None
 |          Maximum size gap to forward or backward fill.
 |      regex : bool or same types as `to_replace`, default False
 |          Whether to interpret `to_replace` and/or `value` as regular
 |          expressions. If this is ``True`` then `to_replace` *must* be a
 |          string. Alternatively, this could be a regular expression or a
 |          list, dict, or array of regular expressions in which case
 |          `to_replace` must be ``None``.
 |      method : {'pad', 'ffill', 'bfill', `None`}
 |          The method to use when for replacement, when `to_replace` is a
 |          scalar, list or tuple and `value` is ``None``.
 |      
 |          .. versionchanged:: 0.23.0
 |              Added to DataFrame.
 |      
 |      See Also
 |      --------
 |      DataFrame.fillna : Fill NA values
 |      DataFrame.where : Replace values based on boolean condition
 |      Series.str.replace : Simple string replacement.
 |      
 |      Returns
 |      -------
 |      DataFrame
 |          Object after replacement.
 |      
 |      Raises
 |      ------
 |      AssertionError
 |          * If `regex` is not a ``bool`` and `to_replace` is not
 |            ``None``.
 |      TypeError
 |          * If `to_replace` is a ``dict`` and `value` is not a ``list``,
 |            ``dict``, ``ndarray``, or ``Series``
 |          * If `to_replace` is ``None`` and `regex` is not compilable
 |            into a regular expression or is a list, dict, ndarray, or
 |            Series.
 |          * When replacing multiple ``bool`` or ``datetime64`` objects and
 |            the arguments to `to_replace` does not match the type of the
 |            value being replaced
 |      ValueError
 |          * If a ``list`` or an ``ndarray`` is passed to `to_replace` and
 |            `value` but they are not the same length.
 |      
 |      Notes
 |      -----
 |      * Regex substitution is performed under the hood with ``re.sub``. The
 |        rules for substitution for ``re.sub`` are the same.
 |      * Regular expressions will only substitute on strings, meaning you
 |        cannot provide, for example, a regular expression matching floating
 |        point numbers and expect the columns in your frame that have a
 |        numeric dtype to be matched. However, if those floating point
 |        numbers *are* strings, then you can do this.
 |      * This method has *a lot* of options. You are encouraged to experiment
 |        and play with this method to gain intuition about how it works.
 |      * When dict is used as the `to_replace` value, it is like
 |        key(s) in the dict are the to_replace part and
 |        value(s) in the dict are the value parameter.
 |      
 |      Examples
 |      --------
 |      
 |      **Scalar `to_replace` and `value`**
 |      
 |      >>> s = pd.Series([0, 1, 2, 3, 4])
 |      >>> s.replace(0, 5)
 |      0    5
 |      1    1
 |      2    2
 |      3    3
 |      4    4
 |      dtype: int64
 |      
 |      >>> df = pd.DataFrame({'A': [0, 1, 2, 3, 4],
 |      ...                    'B': [5, 6, 7, 8, 9],
 |      ...                    'C': ['a', 'b', 'c', 'd', 'e']})
 |      >>> df.replace(0, 5)
 |         A  B  C
 |      0  5  5  a
 |      1  1  6  b
 |      2  2  7  c
 |      3  3  8  d
 |      4  4  9  e
 |      
 |      **List-like `to_replace`**
 |      
 |      >>> df.replace([0, 1, 2, 3], 4)
 |         A  B  C
 |      0  4  5  a
 |      1  4  6  b
 |      2  4  7  c
 |      3  4  8  d
 |      4  4  9  e
 |      
 |      >>> df.replace([0, 1, 2, 3], [4, 3, 2, 1])
 |         A  B  C
 |      0  4  5  a
 |      1  3  6  b
 |      2  2  7  c
 |      3  1  8  d
 |      4  4  9  e
 |      
 |      >>> s.replace([1, 2], method='bfill')
 |      0    0
 |      1    3
 |      2    3
 |      3    3
 |      4    4
 |      dtype: int64
 |      
 |      **dict-like `to_replace`**
 |      
 |      >>> df.replace({0: 10, 1: 100})
 |           A  B  C
 |      0   10  5  a
 |      1  100  6  b
 |      2    2  7  c
 |      3    3  8  d
 |      4    4  9  e
 |      
 |      >>> df.replace({'A': 0, 'B': 5}, 100)
 |           A    B  C
 |      0  100  100  a
 |      1    1    6  b
 |      2    2    7  c
 |      3    3    8  d
 |      4    4    9  e
 |      
 |      >>> df.replace({'A': {0: 100, 4: 400}})
 |           A  B  C
 |      0  100  5  a
 |      1    1  6  b
 |      2    2  7  c
 |      3    3  8  d
 |      4  400  9  e
 |      
 |      **Regular expression `to_replace`**
 |      
 |      >>> df = pd.DataFrame({'A': ['bat', 'foo', 'bait'],
 |      ...                    'B': ['abc', 'bar', 'xyz']})
 |      >>> df.replace(to_replace=r'^ba.$', value='new', regex=True)
 |            A    B
 |      0   new  abc
 |      1   foo  new
 |      2  bait  xyz
 |      
 |      >>> df.replace({'A': r'^ba.$'}, {'A': 'new'}, regex=True)
 |            A    B
 |      0   new  abc
 |      1   foo  bar
 |      2  bait  xyz
 |      
 |      >>> df.replace(regex=r'^ba.$', value='new')
 |            A    B
 |      0   new  abc
 |      1   foo  new
 |      2  bait  xyz
 |      
 |      >>> df.replace(regex={r'^ba.$':'new', 'foo':'xyz'})
 |            A    B
 |      0   new  abc
 |      1   xyz  new
 |      2  bait  xyz
 |      
 |      >>> df.replace(regex=[r'^ba.$', 'foo'], value='new')
 |            A    B
 |      0   new  abc
 |      1   new  new
 |      2  bait  xyz
 |      
 |      Note that when replacing multiple ``bool`` or ``datetime64`` objects,
 |      the data types in the `to_replace` parameter must match the data
 |      type of the value being replaced:
 |      
 |      >>> df = pd.DataFrame({'A': [True, False, True],
 |      ...                    'B': [False, True, False]})
 |      >>> df.replace({'a string': 'new value', True: False})  # raises
 |      Traceback (most recent call last):
 |          ...
 |      TypeError: Cannot compare types 'ndarray(dtype=bool)' and 'str'
 |      
 |      This raises a ``TypeError`` because one of the ``dict`` keys is not of
 |      the correct type for replacement.
 |      
 |      Compare the behavior of ``s.replace({'a': None})`` and
 |      ``s.replace('a', None)`` to understand the pecularities
 |      of the `to_replace` parameter:
 |      
 |      >>> s = pd.Series([10, 'a', 'a', 'b', 'a'])
 |      
 |      When one uses a dict as the `to_replace` value, it is like the
 |      value(s) in the dict are equal to the `value` parameter.
 |      ``s.replace({'a': None})`` is equivalent to
 |      ``s.replace(to_replace={'a': None}, value=None, method=None)``:
 |      
 |      >>> s.replace({'a': None})
 |      0      10
 |      1    None
 |      2    None
 |      3       b
 |      4    None
 |      dtype: object
 |      
 |      When ``value=None`` and `to_replace` is a scalar, list or
 |      tuple, `replace` uses the method parameter (default 'pad') to do the
 |      replacement. So this is why the 'a' values are being replaced by 10
 |      in rows 1 and 2 and 'b' in row 4 in this case.
 |      The command ``s.replace('a', None)`` is actually equivalent to
 |      ``s.replace(to_replace='a', value=None, method='pad')``:
 |      
 |      >>> s.replace('a', None)
 |      0    10
 |      1    10
 |      2    10
 |      3     b
 |      4     b
 |      dtype: object
 |  
 |  reset_index(self, level=None, drop=False, inplace=False, col_level=0, col_fill='')
 |      For DataFrame with multi-level index, return new DataFrame with
 |      labeling information in the columns under the index names, defaulting
 |      to 'level_0', 'level_1', etc. if any are None. For a standard index,
 |      the index name will be used (if set), otherwise a default 'index' or
 |      'level_0' (if 'index' is already taken) will be used.
 |      
 |      Parameters
 |      ----------
 |      level : int, str, tuple, or list, default None
 |          Only remove the given levels from the index. Removes all levels by
 |          default
 |      drop : boolean, default False
 |          Do not try to insert index into dataframe columns. This resets
 |          the index to the default integer index.
 |      inplace : boolean, default False
 |          Modify the DataFrame in place (do not create a new object)
 |      col_level : int or str, default 0
 |          If the columns have multiple levels, determines which level the
 |          labels are inserted into. By default it is inserted into the first
 |          level.
 |      col_fill : object, default ''
 |          If the columns have multiple levels, determines how the other
 |          levels are named. If None then the index name is repeated.
 |      
 |      Returns
 |      -------
 |      resetted : DataFrame
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame([('bird',    389.0),
 |      ...                    ('bird',     24.0),
 |      ...                    ('mammal',   80.5),
 |      ...                    ('mammal', np.nan)],
 |      ...                   index=['falcon', 'parrot', 'lion', 'monkey'],
 |      ...                   columns=('class', 'max_speed'))
 |      >>> df
 |               class  max_speed
 |      falcon    bird      389.0
 |      parrot    bird       24.0
 |      lion    mammal       80.5
 |      monkey  mammal        NaN
 |      
 |      When we reset the index, the old index is added as a column, and a
 |      new sequential index is used:
 |      
 |      >>> df.reset_index()
 |          index   class  max_speed
 |      0  falcon    bird      389.0
 |      1  parrot    bird       24.0
 |      2    lion  mammal       80.5
 |      3  monkey  mammal        NaN
 |      
 |      We can use the `drop` parameter to avoid the old index being added as
 |      a column:
 |      
 |      >>> df.reset_index(drop=True)
 |          class  max_speed
 |      0    bird      389.0
 |      1    bird       24.0
 |      2  mammal       80.5
 |      3  mammal        NaN
 |      
 |      You can also use `reset_index` with `MultiIndex`.
 |      
 |      >>> index = pd.MultiIndex.from_tuples([('bird', 'falcon'),
 |      ...                                    ('bird', 'parrot'),
 |      ...                                    ('mammal', 'lion'),
 |      ...                                    ('mammal', 'monkey')],
 |      ...                                   names=['class', 'name'])
 |      >>> columns = pd.MultiIndex.from_tuples([('speed', 'max'),
 |      ...                                      ('species', 'type')])
 |      >>> df = pd.DataFrame([(389.0, 'fly'),
 |      ...                    ( 24.0, 'fly'),
 |      ...                    ( 80.5, 'run'),
 |      ...                    (np.nan, 'jump')],
 |      ...                   index=index,
 |      ...                   columns=columns)
 |      >>> df
 |                     speed species
 |                       max    type
 |      class  name
 |      bird   falcon  389.0     fly
 |             parrot   24.0     fly
 |      mammal lion     80.5     run
 |             monkey    NaN    jump
 |      
 |      If the index has multiple levels, we can reset a subset of them:
 |      
 |      >>> df.reset_index(level='class')
 |               class  speed species
 |                        max    type
 |      name
 |      falcon    bird  389.0     fly
 |      parrot    bird   24.0     fly
 |      lion    mammal   80.5     run
 |      monkey  mammal    NaN    jump
 |      
 |      If we are not dropping the index, by default, it is placed in the top
 |      level. We can place it in another level:
 |      
 |      >>> df.reset_index(level='class', col_level=1)
 |                      speed species
 |               class    max    type
 |      name
 |      falcon    bird  389.0     fly
 |      parrot    bird   24.0     fly
 |      lion    mammal   80.5     run
 |      monkey  mammal    NaN    jump
 |      
 |      When the index is inserted under another level, we can specify under
 |      which one with the parameter `col_fill`:
 |      
 |      >>> df.reset_index(level='class', col_level=1, col_fill='species')
 |                    species  speed species
 |                      class    max    type
 |      name
 |      falcon           bird  389.0     fly
 |      parrot           bird   24.0     fly
 |      lion           mammal   80.5     run
 |      monkey         mammal    NaN    jump
 |      
 |      If we specify a nonexistent level for `col_fill`, it is created:
 |      
 |      >>> df.reset_index(level='class', col_level=1, col_fill='genus')
 |                      genus  speed species
 |                      class    max    type
 |      name
 |      falcon           bird  389.0     fly
 |      parrot           bird   24.0     fly
 |      lion           mammal   80.5     run
 |      monkey         mammal    NaN    jump
 |  
 |  rfloordiv(self, other, axis='columns', level=None, fill_value=None)
 |      Integer division of dataframe and other, element-wise (binary operator `rfloordiv`).
 |      
 |      Equivalent to ``other // dataframe``, but with support to substitute a fill_value for
 |      missing data in one of the inputs.
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      None
 |      
 |      See also
 |      --------
 |      DataFrame.floordiv
 |  
 |  rmod(self, other, axis='columns', level=None, fill_value=None)
 |      Modulo of dataframe and other, element-wise (binary operator `rmod`).
 |      
 |      Equivalent to ``other % dataframe``, but with support to substitute a fill_value for
 |      missing data in one of the inputs.
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      None
 |      
 |      See also
 |      --------
 |      DataFrame.mod
 |  
 |  rmul(self, other, axis='columns', level=None, fill_value=None)
 |      Multiplication of dataframe and other, element-wise (binary operator `rmul`).
 |      
 |      Equivalent to ``other * dataframe``, but with support to substitute a fill_value for
 |      missing data in one of the inputs.
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      None
 |      
 |      See also
 |      --------
 |      DataFrame.mul
 |  
 |  rolling(self, window, min_periods=None, center=False, win_type=None, on=None, axis=0, closed=None)
 |      Provides rolling window calculations.
 |      
 |      .. versionadded:: 0.18.0
 |      
 |      Parameters
 |      ----------
 |      window : int, or offset
 |          Size of the moving window. This is the number of observations used for
 |          calculating the statistic. Each window will be a fixed size.
 |      
 |          If its an offset then this will be the time period of each window. Each
 |          window will be a variable sized based on the observations included in
 |          the time-period. This is only valid for datetimelike indexes. This is
 |          new in 0.19.0
 |      min_periods : int, default None
 |          Minimum number of observations in window required to have a value
 |          (otherwise result is NA). For a window that is specified by an offset,
 |          this will default to 1.
 |      center : boolean, default False
 |          Set the labels at the center of the window.
 |      win_type : string, default None
 |          Provide a window type. If ``None``, all points are evenly weighted.
 |          See the notes below for further information.
 |      on : string, optional
 |          For a DataFrame, column on which to calculate
 |          the rolling window, rather than the index
 |      closed : string, default None
 |          Make the interval closed on the 'right', 'left', 'both' or
 |          'neither' endpoints.
 |          For offset-based windows, it defaults to 'right'.
 |          For fixed windows, defaults to 'both'. Remaining cases not implemented
 |          for fixed windows.
 |      
 |          .. versionadded:: 0.20.0
 |      
 |      axis : int or string, default 0
 |      
 |      Returns
 |      -------
 |      a Window or Rolling sub-classed for the particular operation
 |      
 |      Examples
 |      --------
 |      
 |      >>> df = pd.DataFrame({'B': [0, 1, 2, np.nan, 4]})
 |      >>> df
 |           B
 |      0  0.0
 |      1  1.0
 |      2  2.0
 |      3  NaN
 |      4  4.0
 |      
 |      Rolling sum with a window length of 2, using the 'triang'
 |      window type.
 |      
 |      >>> df.rolling(2, win_type='triang').sum()
 |           B
 |      0  NaN
 |      1  1.0
 |      2  2.5
 |      3  NaN
 |      4  NaN
 |      
 |      Rolling sum with a window length of 2, min_periods defaults
 |      to the window length.
 |      
 |      >>> df.rolling(2).sum()
 |           B
 |      0  NaN
 |      1  1.0
 |      2  3.0
 |      3  NaN
 |      4  NaN
 |      
 |      Same as above, but explicitly set the min_periods
 |      
 |      >>> df.rolling(2, min_periods=1).sum()
 |           B
 |      0  0.0
 |      1  1.0
 |      2  3.0
 |      3  2.0
 |      4  4.0
 |      
 |      A ragged (meaning not-a-regular frequency), time-indexed DataFrame
 |      
 |      >>> df = pd.DataFrame({'B': [0, 1, 2, np.nan, 4]},
 |      ...                   index = [pd.Timestamp('20130101 09:00:00'),
 |      ...                            pd.Timestamp('20130101 09:00:02'),
 |      ...                            pd.Timestamp('20130101 09:00:03'),
 |      ...                            pd.Timestamp('20130101 09:00:05'),
 |      ...                            pd.Timestamp('20130101 09:00:06')])
 |      
 |      >>> df
 |                             B
 |      2013-01-01 09:00:00  0.0
 |      2013-01-01 09:00:02  1.0
 |      2013-01-01 09:00:03  2.0
 |      2013-01-01 09:00:05  NaN
 |      2013-01-01 09:00:06  4.0
 |      
 |      
 |      Contrasting to an integer rolling window, this will roll a variable
 |      length window corresponding to the time period.
 |      The default for min_periods is 1.
 |      
 |      >>> df.rolling('2s').sum()
 |                             B
 |      2013-01-01 09:00:00  0.0
 |      2013-01-01 09:00:02  1.0
 |      2013-01-01 09:00:03  3.0
 |      2013-01-01 09:00:05  NaN
 |      2013-01-01 09:00:06  4.0
 |      
 |      Notes
 |      -----
 |      By default, the result is set to the right edge of the window. This can be
 |      changed to the center of the window by setting ``center=True``.
 |      
 |      To learn more about the offsets & frequency strings, please see `this link
 |      <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__.
 |      
 |      The recognized win_types are:
 |      
 |      * ``boxcar``
 |      * ``triang``
 |      * ``blackman``
 |      * ``hamming``
 |      * ``bartlett``
 |      * ``parzen``
 |      * ``bohman``
 |      * ``blackmanharris``
 |      * ``nuttall``
 |      * ``barthann``
 |      * ``kaiser`` (needs beta)
 |      * ``gaussian`` (needs std)
 |      * ``general_gaussian`` (needs power, width)
 |      * ``slepian`` (needs width).
 |      
 |      If ``win_type=None`` all points are evenly weighted. To learn more about
 |      different window types see `scipy.signal window functions
 |      <https://docs.scipy.org/doc/scipy/reference/signal.html#window-functions>`__.
 |      
 |      See Also
 |      --------
 |      expanding : Provides expanding transformations.
 |      ewm : Provides exponential weighted functions
 |  
 |  round(self, decimals=0, *args, **kwargs)
 |      Round a DataFrame to a variable number of decimal places.
 |      
 |      Parameters
 |      ----------
 |      decimals : int, dict, Series
 |          Number of decimal places to round each column to. If an int is
 |          given, round each column to the same number of places.
 |          Otherwise dict and Series round to variable numbers of places.
 |          Column names should be in the keys if `decimals` is a
 |          dict-like, or in the index if `decimals` is a Series. Any
 |          columns not included in `decimals` will be left as is. Elements
 |          of `decimals` which are not columns of the input will be
 |          ignored.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame(np.random.random([3, 3]),
 |      ...     columns=['A', 'B', 'C'], index=['first', 'second', 'third'])
 |      >>> df
 |                     A         B         C
 |      first   0.028208  0.992815  0.173891
 |      second  0.038683  0.645646  0.577595
 |      third   0.877076  0.149370  0.491027
 |      >>> df.round(2)
 |                 A     B     C
 |      first   0.03  0.99  0.17
 |      second  0.04  0.65  0.58
 |      third   0.88  0.15  0.49
 |      >>> df.round({'A': 1, 'C': 2})
 |                A         B     C
 |      first   0.0  0.992815  0.17
 |      second  0.0  0.645646  0.58
 |      third   0.9  0.149370  0.49
 |      >>> decimals = pd.Series([1, 0, 2], index=['A', 'B', 'C'])
 |      >>> df.round(decimals)
 |                A  B     C
 |      first   0.0  1  0.17
 |      second  0.0  1  0.58
 |      third   0.9  0  0.49
 |      
 |      Returns
 |      -------
 |      DataFrame object
 |      
 |      See Also
 |      --------
 |      numpy.around
 |      Series.round
 |  
 |  rpow(self, other, axis='columns', level=None, fill_value=None)
 |      Exponential power of dataframe and other, element-wise (binary operator `rpow`).
 |      
 |      Equivalent to ``other ** dataframe``, but with support to substitute a fill_value for
 |      missing data in one of the inputs.
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      None
 |      
 |      See also
 |      --------
 |      DataFrame.pow
 |  
 |  rsub(self, other, axis='columns', level=None, fill_value=None)
 |      Subtraction of dataframe and other, element-wise (binary operator `rsub`).
 |      
 |      Equivalent to ``other - dataframe``, but with support to substitute a fill_value for
 |      missing data in one of the inputs.
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      
 |      >>> a = pd.DataFrame([2, 1, 1, np.nan], index=['a', 'b', 'c', 'd'],
 |      ...                  columns=['one'])
 |      >>> a
 |         one
 |      a  2.0
 |      b  1.0
 |      c  1.0
 |      d  NaN
 |      >>> b = pd.DataFrame(dict(one=[1, np.nan, 1, np.nan],
 |      ...                       two=[3, 2, np.nan, 2]),
 |      ...                  index=['a', 'b', 'd', 'e'])
 |      >>> b
 |         one  two
 |      a  1.0  3.0
 |      b  NaN  2.0
 |      d  1.0  NaN
 |      e  NaN  2.0
 |      >>> a.sub(b, fill_value=0)
 |         one  two
 |      a  1.0  -3.0
 |      b  1.0  -2.0
 |      c  1.0  NaN
 |      d  -1.0  NaN
 |      e  NaN  -2.0
 |      
 |      
 |      See also
 |      --------
 |      DataFrame.sub
 |  
 |  rtruediv(self, other, axis='columns', level=None, fill_value=None)
 |      Floating division of dataframe and other, element-wise (binary operator `rtruediv`).
 |      
 |      Equivalent to ``other / dataframe``, but with support to substitute a fill_value for
 |      missing data in one of the inputs.
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      None
 |      
 |      See also
 |      --------
 |      DataFrame.truediv
 |  
 |  select_dtypes(self, include=None, exclude=None)
 |      Return a subset of the DataFrame's columns based on the column dtypes.
 |      
 |      Parameters
 |      ----------
 |      include, exclude : scalar or list-like
 |          A selection of dtypes or strings to be included/excluded. At least
 |          one of these parameters must be supplied.
 |      
 |      Raises
 |      ------
 |      ValueError
 |          * If both of ``include`` and ``exclude`` are empty
 |          * If ``include`` and ``exclude`` have overlapping elements
 |          * If any kind of string dtype is passed in.
 |      
 |      Returns
 |      -------
 |      subset : DataFrame
 |          The subset of the frame including the dtypes in ``include`` and
 |          excluding the dtypes in ``exclude``.
 |      
 |      Notes
 |      -----
 |      * To select all *numeric* types, use ``np.number`` or ``'number'``
 |      * To select strings you must use the ``object`` dtype, but note that
 |        this will return *all* object dtype columns
 |      * See the `numpy dtype hierarchy
 |        <http://docs.scipy.org/doc/numpy/reference/arrays.scalars.html>`__
 |      * To select datetimes, use ``np.datetime64``, ``'datetime'`` or
 |        ``'datetime64'``
 |      * To select timedeltas, use ``np.timedelta64``, ``'timedelta'`` or
 |        ``'timedelta64'``
 |      * To select Pandas categorical dtypes, use ``'category'``
 |      * To select Pandas datetimetz dtypes, use ``'datetimetz'`` (new in
 |        0.20.0) or ``'datetime64[ns, tz]'``
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'a': [1, 2] * 3,
 |      ...                    'b': [True, False] * 3,
 |      ...                    'c': [1.0, 2.0] * 3})
 |      >>> df
 |              a      b  c
 |      0       1   True  1.0
 |      1       2  False  2.0
 |      2       1   True  1.0
 |      3       2  False  2.0
 |      4       1   True  1.0
 |      5       2  False  2.0
 |      
 |      >>> df.select_dtypes(include='bool')
 |         b
 |      0  True
 |      1  False
 |      2  True
 |      3  False
 |      4  True
 |      5  False
 |      
 |      >>> df.select_dtypes(include=['float64'])
 |         c
 |      0  1.0
 |      1  2.0
 |      2  1.0
 |      3  2.0
 |      4  1.0
 |      5  2.0
 |      
 |      >>> df.select_dtypes(exclude=['int'])
 |             b    c
 |      0   True  1.0
 |      1  False  2.0
 |      2   True  1.0
 |      3  False  2.0
 |      4   True  1.0
 |      5  False  2.0
 |  
 |  sem(self, axis=None, skipna=None, level=None, ddof=1, numeric_only=None, **kwargs)
 |      Return unbiased standard error of the mean over requested axis.
 |      
 |      Normalized by N-1 by default. This can be changed using the ddof argument
 |      
 |      Parameters
 |      ----------
 |      axis : {index (0), columns (1)}
 |      skipna : boolean, default True
 |          Exclude NA/null values. If an entire row/column is NA, the result
 |          will be NA
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series
 |      ddof : int, default 1
 |          Delta Degrees of Freedom. The divisor used in calculations is N - ddof,
 |          where N represents the number of elements.
 |      numeric_only : boolean, default None
 |          Include only float, int, boolean columns. If None, will attempt to use
 |          everything, then use only numeric data. Not implemented for Series.
 |      
 |      Returns
 |      -------
 |      sem : Series or DataFrame (if level specified)
 |  
 |  set_index(self, keys, drop=True, append=False, inplace=False, verify_integrity=False)
 |      Set the DataFrame index (row labels) using one or more existing
 |      columns. By default yields a new object.
 |      
 |      Parameters
 |      ----------
 |      keys : column label or list of column labels / arrays
 |      drop : boolean, default True
 |          Delete columns to be used as the new index
 |      append : boolean, default False
 |          Whether to append columns to existing index
 |      inplace : boolean, default False
 |          Modify the DataFrame in place (do not create a new object)
 |      verify_integrity : boolean, default False
 |          Check the new index for duplicates. Otherwise defer the check until
 |          necessary. Setting to False will improve the performance of this
 |          method
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'month': [1, 4, 7, 10],
 |      ...                    'year': [2012, 2014, 2013, 2014],
 |      ...                    'sale':[55, 40, 84, 31]})
 |         month  sale  year
 |      0  1      55    2012
 |      1  4      40    2014
 |      2  7      84    2013
 |      3  10     31    2014
 |      
 |      Set the index to become the 'month' column:
 |      
 |      >>> df.set_index('month')
 |             sale  year
 |      month
 |      1      55    2012
 |      4      40    2014
 |      7      84    2013
 |      10     31    2014
 |      
 |      Create a multi-index using columns 'year' and 'month':
 |      
 |      >>> df.set_index(['year', 'month'])
 |                  sale
 |      year  month
 |      2012  1     55
 |      2014  4     40
 |      2013  7     84
 |      2014  10    31
 |      
 |      Create a multi-index using a set of values and a column:
 |      
 |      >>> df.set_index([[1, 2, 3, 4], 'year'])
 |               month  sale
 |         year
 |      1  2012  1      55
 |      2  2014  4      40
 |      3  2013  7      84
 |      4  2014  10     31
 |      
 |      Returns
 |      -------
 |      dataframe : DataFrame
 |  
 |  set_value(self, index, col, value, takeable=False)
 |      Put single value at passed column and index
 |      
 |      .. deprecated:: 0.21.0
 |          Use .at[] or .iat[] accessors instead.
 |      
 |      Parameters
 |      ----------
 |      index : row label
 |      col : column label
 |      value : scalar value
 |      takeable : interpret the index/col as indexers, default False
 |      
 |      Returns
 |      -------
 |      frame : DataFrame
 |          If label pair is contained, will be reference to calling DataFrame,
 |          otherwise a new object
 |  
 |  shift(self, periods=1, freq=None, axis=0)
 |      Shift index by desired number of periods with an optional time freq
 |      
 |      Parameters
 |      ----------
 |      periods : int
 |          Number of periods to move, can be positive or negative
 |      freq : DateOffset, timedelta, or time rule string, optional
 |          Increment to use from the tseries module or time rule (e.g. 'EOM').
 |          See Notes.
 |      axis : {0 or 'index', 1 or 'columns'}
 |      
 |      Notes
 |      -----
 |      If freq is specified then the index values are shifted but the data
 |      is not realigned. That is, use freq if you would like to extend the
 |      index when shifting and preserve the original data.
 |      
 |      Returns
 |      -------
 |      shifted : DataFrame
 |  
 |  skew(self, axis=None, skipna=None, level=None, numeric_only=None, **kwargs)
 |      Return unbiased skew over requested axis
 |      Normalized by N-1
 |      
 |      Parameters
 |      ----------
 |      axis : {index (0), columns (1)}
 |      skipna : boolean, default True
 |          Exclude NA/null values when computing the result.
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series
 |      numeric_only : boolean, default None
 |          Include only float, int, boolean columns. If None, will attempt to use
 |          everything, then use only numeric data. Not implemented for Series.
 |      
 |      Returns
 |      -------
 |      skew : Series or DataFrame (if level specified)
 |  
 |  sort_index(self, axis=0, level=None, ascending=True, inplace=False, kind='quicksort', na_position='last', sort_remaining=True, by=None)
 |      Sort object by labels (along an axis)
 |      
 |      Parameters
 |      ----------
 |      axis : index, columns to direct sorting
 |      level : int or level name or list of ints or list of level names
 |          if not None, sort on values in specified index level(s)
 |      ascending : boolean, default True
 |          Sort ascending vs. descending
 |      inplace : bool, default False
 |          if True, perform operation in-place
 |      kind : {'quicksort', 'mergesort', 'heapsort'}, default 'quicksort'
 |           Choice of sorting algorithm. See also ndarray.np.sort for more
 |           information.  `mergesort` is the only stable algorithm. For
 |           DataFrames, this option is only applied when sorting on a single
 |           column or label.
 |      na_position : {'first', 'last'}, default 'last'
 |           `first` puts NaNs at the beginning, `last` puts NaNs at the end.
 |           Not implemented for MultiIndex.
 |      sort_remaining : bool, default True
 |          if true and sorting by level and index is multilevel, sort by other
 |          levels too (in order) after sorting by specified level
 |      
 |      Returns
 |      -------
 |      sorted_obj : DataFrame
 |  
 |  sort_values(self, by, axis=0, ascending=True, inplace=False, kind='quicksort', na_position='last')
 |      Sort by the values along either axis
 |      
 |      Parameters
 |      ----------
 |      by : str or list of str
 |          Name or list of names to sort by.
 |      
 |          - if `axis` is 0 or `'index'` then `by` may contain index
 |            levels and/or column labels
 |          - if `axis` is 1 or `'columns'` then `by` may contain column
 |            levels and/or index labels
 |      
 |          .. versionchanged:: 0.23.0
 |             Allow specifying index or column level names.
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |           Axis to be sorted
 |      ascending : bool or list of bool, default True
 |           Sort ascending vs. descending. Specify list for multiple sort
 |           orders.  If this is a list of bools, must match the length of
 |           the by.
 |      inplace : bool, default False
 |           if True, perform operation in-place
 |      kind : {'quicksort', 'mergesort', 'heapsort'}, default 'quicksort'
 |           Choice of sorting algorithm. See also ndarray.np.sort for more
 |           information.  `mergesort` is the only stable algorithm. For
 |           DataFrames, this option is only applied when sorting on a single
 |           column or label.
 |      na_position : {'first', 'last'}, default 'last'
 |           `first` puts NaNs at the beginning, `last` puts NaNs at the end
 |      
 |      Returns
 |      -------
 |      sorted_obj : DataFrame
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({
 |      ...     'col1' : ['A', 'A', 'B', np.nan, 'D', 'C'],
 |      ...     'col2' : [2, 1, 9, 8, 7, 4],
 |      ...     'col3': [0, 1, 9, 4, 2, 3],
 |      ... })
 |      >>> df
 |          col1 col2 col3
 |      0   A    2    0
 |      1   A    1    1
 |      2   B    9    9
 |      3   NaN  8    4
 |      4   D    7    2
 |      5   C    4    3
 |      
 |      Sort by col1
 |      
 |      >>> df.sort_values(by=['col1'])
 |          col1 col2 col3
 |      0   A    2    0
 |      1   A    1    1
 |      2   B    9    9
 |      5   C    4    3
 |      4   D    7    2
 |      3   NaN  8    4
 |      
 |      Sort by multiple columns
 |      
 |      >>> df.sort_values(by=['col1', 'col2'])
 |          col1 col2 col3
 |      1   A    1    1
 |      0   A    2    0
 |      2   B    9    9
 |      5   C    4    3
 |      4   D    7    2
 |      3   NaN  8    4
 |      
 |      Sort Descending
 |      
 |      >>> df.sort_values(by='col1', ascending=False)
 |          col1 col2 col3
 |      4   D    7    2
 |      5   C    4    3
 |      2   B    9    9
 |      0   A    2    0
 |      1   A    1    1
 |      3   NaN  8    4
 |      
 |      Putting NAs first
 |      
 |      >>> df.sort_values(by='col1', ascending=False, na_position='first')
 |          col1 col2 col3
 |      3   NaN  8    4
 |      4   D    7    2
 |      5   C    4    3
 |      2   B    9    9
 |      0   A    2    0
 |      1   A    1    1
 |  
 |  sortlevel(self, level=0, axis=0, ascending=True, inplace=False, sort_remaining=True)
 |      Sort multilevel index by chosen axis and primary level. Data will be
 |      lexicographically sorted by the chosen level followed by the other
 |      levels (in order).
 |      
 |      .. deprecated:: 0.20.0
 |          Use :meth:`DataFrame.sort_index`
 |      
 |      
 |      Parameters
 |      ----------
 |      level : int
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |      ascending : boolean, default True
 |      inplace : boolean, default False
 |          Sort the DataFrame without creating a new instance
 |      sort_remaining : boolean, default True
 |          Sort by the other levels too.
 |      
 |      Returns
 |      -------
 |      sorted : DataFrame
 |      
 |      See Also
 |      --------
 |      DataFrame.sort_index(level=...)
 |  
 |  stack(self, level=-1, dropna=True)
 |      Stack the prescribed level(s) from columns to index.
 |      
 |      Return a reshaped DataFrame or Series having a multi-level
 |      index with one or more new inner-most levels compared to the current
 |      DataFrame. The new inner-most levels are created by pivoting the
 |      columns of the current dataframe:
 |      
 |        - if the columns have a single level, the output is a Series;
 |        - if the columns have multiple levels, the new index
 |          level(s) is (are) taken from the prescribed level(s) and
 |          the output is a DataFrame.
 |      
 |      The new index levels are sorted.
 |      
 |      Parameters
 |      ----------
 |      level : int, str, list, default -1
 |          Level(s) to stack from the column axis onto the index
 |          axis, defined as one index or label, or a list of indices
 |          or labels.
 |      dropna : bool, default True
 |          Whether to drop rows in the resulting Frame/Series with
 |          missing values. Stacking a column level onto the index
 |          axis can create combinations of index and column values
 |          that are missing from the original dataframe. See Examples
 |          section.
 |      
 |      Returns
 |      -------
 |      DataFrame or Series
 |          Stacked dataframe or series.
 |      
 |      See Also
 |      --------
 |      DataFrame.unstack : Unstack prescribed level(s) from index axis
 |           onto column axis.
 |      DataFrame.pivot : Reshape dataframe from long format to wide
 |           format.
 |      DataFrame.pivot_table : Create a spreadsheet-style pivot table
 |           as a DataFrame.
 |      
 |      Notes
 |      -----
 |      The function is named by analogy with a collection of books
 |      being re-organised from being side by side on a horizontal
 |      position (the columns of the dataframe) to being stacked
 |      vertically on top of of each other (in the index of the
 |      dataframe).
 |      
 |      Examples
 |      --------
 |      **Single level columns**
 |      
 |      >>> df_single_level_cols = pd.DataFrame([[0, 1], [2, 3]],
 |      ...                                     index=['cat', 'dog'],
 |      ...                                     columns=['weight', 'height'])
 |      
 |      Stacking a dataframe with a single level column axis returns a Series:
 |      
 |      >>> df_single_level_cols
 |           weight height
 |      cat       0      1
 |      dog       2      3
 |      >>> df_single_level_cols.stack()
 |      cat  weight    0
 |           height    1
 |      dog  weight    2
 |           height    3
 |      dtype: int64
 |      
 |      **Multi level columns: simple case**
 |      
 |      >>> multicol1 = pd.MultiIndex.from_tuples([('weight', 'kg'),
 |      ...                                        ('weight', 'pounds')])
 |      >>> df_multi_level_cols1 = pd.DataFrame([[1, 2], [2, 4]],
 |      ...                                     index=['cat', 'dog'],
 |      ...                                     columns=multicol1)
 |      
 |      Stacking a dataframe with a multi-level column axis:
 |      
 |      >>> df_multi_level_cols1
 |           weight
 |               kg    pounds
 |      cat       1        2
 |      dog       2        4
 |      >>> df_multi_level_cols1.stack()
 |                  weight
 |      cat kg           1
 |          pounds       2
 |      dog kg           2
 |          pounds       4
 |      
 |      **Missing values**
 |      
 |      >>> multicol2 = pd.MultiIndex.from_tuples([('weight', 'kg'),
 |      ...                                        ('height', 'm')])
 |      >>> df_multi_level_cols2 = pd.DataFrame([[1.0, 2.0], [3.0, 4.0]],
 |      ...                                     index=['cat', 'dog'],
 |      ...                                     columns=multicol2)
 |      
 |      It is common to have missing values when stacking a dataframe
 |      with multi-level columns, as the stacked dataframe typically
 |      has more values than the original dataframe. Missing values
 |      are filled with NaNs:
 |      
 |      >>> df_multi_level_cols2
 |          weight height
 |              kg      m
 |      cat    1.0    2.0
 |      dog    3.0    4.0
 |      >>> df_multi_level_cols2.stack()
 |              height  weight
 |      cat kg     NaN     1.0
 |          m      2.0     NaN
 |      dog kg     NaN     3.0
 |          m      4.0     NaN
 |      
 |      **Prescribing the level(s) to be stacked**
 |      
 |      The first parameter controls which level or levels are stacked:
 |      
 |      >>> df_multi_level_cols2.stack(0)
 |                   kg    m
 |      cat height  NaN  2.0
 |          weight  1.0  NaN
 |      dog height  NaN  4.0
 |          weight  3.0  NaN
 |      >>> df_multi_level_cols2.stack([0, 1])
 |      cat  height  m     2.0
 |           weight  kg    1.0
 |      dog  height  m     4.0
 |           weight  kg    3.0
 |      dtype: float64
 |      
 |      **Dropping missing values**
 |      
 |      >>> df_multi_level_cols3 = pd.DataFrame([[None, 1.0], [2.0, 3.0]],
 |      ...                                     index=['cat', 'dog'],
 |      ...                                     columns=multicol2)
 |      
 |      Note that rows where all values are missing are dropped by
 |      default but this behaviour can be controlled via the dropna
 |      keyword parameter:
 |      
 |      >>> df_multi_level_cols3
 |          weight height
 |              kg      m
 |      cat    NaN    1.0
 |      dog    2.0    3.0
 |      >>> df_multi_level_cols3.stack(dropna=False)
 |              height  weight
 |      cat kg     NaN     NaN
 |          m      1.0     NaN
 |      dog kg     NaN     2.0
 |          m      3.0     NaN
 |      >>> df_multi_level_cols3.stack(dropna=True)
 |              height  weight
 |      cat m      1.0     NaN
 |      dog kg     NaN     2.0
 |          m      3.0     NaN
 |  
 |  std(self, axis=None, skipna=None, level=None, ddof=1, numeric_only=None, **kwargs)
 |      Return sample standard deviation over requested axis.
 |      
 |      Normalized by N-1 by default. This can be changed using the ddof argument
 |      
 |      Parameters
 |      ----------
 |      axis : {index (0), columns (1)}
 |      skipna : boolean, default True
 |          Exclude NA/null values. If an entire row/column is NA, the result
 |          will be NA
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series
 |      ddof : int, default 1
 |          Delta Degrees of Freedom. The divisor used in calculations is N - ddof,
 |          where N represents the number of elements.
 |      numeric_only : boolean, default None
 |          Include only float, int, boolean columns. If None, will attempt to use
 |          everything, then use only numeric data. Not implemented for Series.
 |      
 |      Returns
 |      -------
 |      std : Series or DataFrame (if level specified)
 |  
 |  sub(self, other, axis='columns', level=None, fill_value=None)
 |      Subtraction of dataframe and other, element-wise (binary operator `sub`).
 |      
 |      Equivalent to ``dataframe - other``, but with support to substitute a fill_value for
 |      missing data in one of the inputs.
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      
 |      >>> a = pd.DataFrame([2, 1, 1, np.nan], index=['a', 'b', 'c', 'd'],
 |      ...                  columns=['one'])
 |      >>> a
 |         one
 |      a  2.0
 |      b  1.0
 |      c  1.0
 |      d  NaN
 |      >>> b = pd.DataFrame(dict(one=[1, np.nan, 1, np.nan],
 |      ...                       two=[3, 2, np.nan, 2]),
 |      ...                  index=['a', 'b', 'd', 'e'])
 |      >>> b
 |         one  two
 |      a  1.0  3.0
 |      b  NaN  2.0
 |      d  1.0  NaN
 |      e  NaN  2.0
 |      >>> a.sub(b, fill_value=0)
 |         one  two
 |      a  1.0  -3.0
 |      b  1.0  -2.0
 |      c  1.0  NaN
 |      d  -1.0  NaN
 |      e  NaN  -2.0
 |      
 |      
 |      See also
 |      --------
 |      DataFrame.rsub
 |  
 |  subtract = sub(self, other, axis='columns', level=None, fill_value=None)
 |  
 |  sum(self, axis=None, skipna=None, level=None, numeric_only=None, min_count=0, **kwargs)
 |      Return the sum of the values for the requested axis
 |      
 |      Parameters
 |      ----------
 |      axis : {index (0), columns (1)}
 |      skipna : boolean, default True
 |          Exclude NA/null values when computing the result.
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series
 |      numeric_only : boolean, default None
 |          Include only float, int, boolean columns. If None, will attempt to use
 |          everything, then use only numeric data. Not implemented for Series.
 |      min_count : int, default 0
 |          The required number of valid values to perform the operation. If fewer than
 |          ``min_count`` non-NA values are present the result will be NA.
 |      
 |          .. versionadded :: 0.22.0
 |      
 |             Added with the default being 0. This means the sum of an all-NA
 |             or empty Series is 0, and the product of an all-NA or empty
 |             Series is 1.
 |      
 |      Returns
 |      -------
 |      sum : Series or DataFrame (if level specified)
 |      
 |      Examples
 |      --------
 |      By default, the sum of an empty or all-NA Series is ``0``.
 |      
 |      >>> pd.Series([]).sum()  # min_count=0 is the default
 |      0.0
 |      
 |      This can be controlled with the ``min_count`` parameter. For example, if
 |      you'd like the sum of an empty series to be NaN, pass ``min_count=1``.
 |      
 |      >>> pd.Series([]).sum(min_count=1)
 |      nan
 |      
 |      Thanks to the ``skipna`` parameter, ``min_count`` handles all-NA and
 |      empty series identically.
 |      
 |      >>> pd.Series([np.nan]).sum()
 |      0.0
 |      
 |      >>> pd.Series([np.nan]).sum(min_count=1)
 |      nan
 |  
 |  swaplevel(self, i=-2, j=-1, axis=0)
 |      Swap levels i and j in a MultiIndex on a particular axis
 |      
 |      Parameters
 |      ----------
 |      i, j : int, string (can be mixed)
 |          Level of index to be swapped. Can pass level name as string.
 |      
 |      Returns
 |      -------
 |      swapped : type of caller (new object)
 |      
 |      .. versionchanged:: 0.18.1
 |      
 |         The indexes ``i`` and ``j`` are now optional, and default to
 |         the two innermost levels of the index.
 |  
 |  to_csv(self, path_or_buf=None, sep=',', na_rep='', float_format=None, columns=None, header=True, index=True, index_label=None, mode='w', encoding=None, compression=None, quoting=None, quotechar='"', line_terminator='\n', chunksize=None, tupleize_cols=None, date_format=None, doublequote=True, escapechar=None, decimal='.')
 |      Write DataFrame to a comma-separated values (csv) file
 |      
 |      Parameters
 |      ----------
 |      path_or_buf : string or file handle, default None
 |          File path or object, if None is provided the result is returned as
 |          a string.
 |      sep : character, default ','
 |          Field delimiter for the output file.
 |      na_rep : string, default ''
 |          Missing data representation
 |      float_format : string, default None
 |          Format string for floating point numbers
 |      columns : sequence, optional
 |          Columns to write
 |      header : boolean or list of string, default True
 |          Write out the column names. If a list of strings is given it is
 |          assumed to be aliases for the column names
 |      index : boolean, default True
 |          Write row names (index)
 |      index_label : string or sequence, or False, default None
 |          Column label for index column(s) if desired. If None is given, and
 |          `header` and `index` are True, then the index names are used. A
 |          sequence should be given if the DataFrame uses MultiIndex.  If
 |          False do not print fields for index names. Use index_label=False
 |          for easier importing in R
 |      mode : str
 |          Python write mode, default 'w'
 |      encoding : string, optional
 |          A string representing the encoding to use in the output file,
 |          defaults to 'ascii' on Python 2 and 'utf-8' on Python 3.
 |      compression : string, optional
 |          A string representing the compression to use in the output file.
 |          Allowed values are 'gzip', 'bz2', 'zip', 'xz'. This input is only
 |          used when the first argument is a filename.
 |      line_terminator : string, default ``'\n'``
 |          The newline character or character sequence to use in the output
 |          file
 |      quoting : optional constant from csv module
 |          defaults to csv.QUOTE_MINIMAL. If you have set a `float_format`
 |          then floats are converted to strings and thus csv.QUOTE_NONNUMERIC
 |          will treat them as non-numeric
 |      quotechar : string (length 1), default '\"'
 |          character used to quote fields
 |      doublequote : boolean, default True
 |          Control quoting of `quotechar` inside a field
 |      escapechar : string (length 1), default None
 |          character used to escape `sep` and `quotechar` when appropriate
 |      chunksize : int or None
 |          rows to write at a time
 |      tupleize_cols : boolean, default False
 |          .. deprecated:: 0.21.0
 |             This argument will be removed and will always write each row
 |             of the multi-index as a separate row in the CSV file.
 |      
 |          Write MultiIndex columns as a list of tuples (if True) or in
 |          the new, expanded format, where each MultiIndex column is a row
 |          in the CSV (if False).
 |      date_format : string, default None
 |          Format string for datetime objects
 |      decimal: string, default '.'
 |          Character recognized as decimal separator. E.g. use ',' for
 |          European data
 |  
 |  to_dict(self, orient='dict', into=<class 'dict'>)
 |      Convert the DataFrame to a dictionary.
 |      
 |      The type of the key-value pairs can be customized with the parameters
 |      (see below).
 |      
 |      Parameters
 |      ----------
 |      orient : str {'dict', 'list', 'series', 'split', 'records', 'index'}
 |          Determines the type of the values of the dictionary.
 |      
 |          - 'dict' (default) : dict like {column -> {index -> value}}
 |          - 'list' : dict like {column -> [values]}
 |          - 'series' : dict like {column -> Series(values)}
 |          - 'split' : dict like
 |            {'index' -> [index], 'columns' -> [columns], 'data' -> [values]}
 |          - 'records' : list like
 |            [{column -> value}, ... , {column -> value}]
 |          - 'index' : dict like {index -> {column -> value}}
 |      
 |          Abbreviations are allowed. `s` indicates `series` and `sp`
 |          indicates `split`.
 |      
 |      into : class, default dict
 |          The collections.Mapping subclass used for all Mappings
 |          in the return value.  Can be the actual class or an empty
 |          instance of the mapping type you want.  If you want a
 |          collections.defaultdict, you must pass it initialized.
 |      
 |          .. versionadded:: 0.21.0
 |      
 |      Returns
 |      -------
 |      result : collections.Mapping like {column -> {index -> value}}
 |      
 |      See Also
 |      --------
 |      DataFrame.from_dict: create a DataFrame from a dictionary
 |      DataFrame.to_json: convert a DataFrame to JSON format
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'col1': [1, 2],
 |      ...                    'col2': [0.5, 0.75]},
 |      ...                   index=['a', 'b'])
 |      >>> df
 |         col1  col2
 |      a     1   0.50
 |      b     2   0.75
 |      >>> df.to_dict()
 |      {'col1': {'a': 1, 'b': 2}, 'col2': {'a': 0.5, 'b': 0.75}}
 |      
 |      You can specify the return orientation.
 |      
 |      >>> df.to_dict('series')
 |      {'col1': a    1
 |               b    2
 |               Name: col1, dtype: int64,
 |       'col2': a    0.50
 |               b    0.75
 |               Name: col2, dtype: float64}
 |      
 |      >>> df.to_dict('split')
 |      {'index': ['a', 'b'], 'columns': ['col1', 'col2'],
 |       'data': [[1.0, 0.5], [2.0, 0.75]]}
 |      
 |      >>> df.to_dict('records')
 |      [{'col1': 1.0, 'col2': 0.5}, {'col1': 2.0, 'col2': 0.75}]
 |      
 |      >>> df.to_dict('index')
 |      {'a': {'col1': 1.0, 'col2': 0.5}, 'b': {'col1': 2.0, 'col2': 0.75}}
 |      
 |      You can also specify the mapping type.
 |      
 |      >>> from collections import OrderedDict, defaultdict
 |      >>> df.to_dict(into=OrderedDict)
 |      OrderedDict([('col1', OrderedDict([('a', 1), ('b', 2)])),
 |                   ('col2', OrderedDict([('a', 0.5), ('b', 0.75)]))])
 |      
 |      If you want a `defaultdict`, you need to initialize it:
 |      
 |      >>> dd = defaultdict(list)
 |      >>> df.to_dict('records', into=dd)
 |      [defaultdict(<class 'list'>, {'col1': 1.0, 'col2': 0.5}),
 |       defaultdict(<class 'list'>, {'col1': 2.0, 'col2': 0.75})]
 |  
 |  to_excel(self, excel_writer, sheet_name='Sheet1', na_rep='', float_format=None, columns=None, header=True, index=True, index_label=None, startrow=0, startcol=0, engine=None, merge_cells=True, encoding=None, inf_rep='inf', verbose=True, freeze_panes=None)
 |      Write DataFrame to an excel sheet
 |      
 |      
 |      Parameters
 |      ----------
 |      excel_writer : string or ExcelWriter object
 |          File path or existing ExcelWriter
 |      sheet_name : string, default 'Sheet1'
 |          Name of sheet which will contain DataFrame
 |      na_rep : string, default ''
 |          Missing data representation
 |      float_format : string, default None
 |          Format string for floating point numbers
 |      columns : sequence, optional
 |          Columns to write
 |      header : boolean or list of string, default True
 |          Write out the column names. If a list of strings is given it is
 |          assumed to be aliases for the column names
 |      index : boolean, default True
 |          Write row names (index)
 |      index_label : string or sequence, default None
 |          Column label for index column(s) if desired. If None is given, and
 |          `header` and `index` are True, then the index names are used. A
 |          sequence should be given if the DataFrame uses MultiIndex.
 |      startrow :
 |          upper left cell row to dump data frame
 |      startcol :
 |          upper left cell column to dump data frame
 |      engine : string, default None
 |          write engine to use - you can also set this via the options
 |          ``io.excel.xlsx.writer``, ``io.excel.xls.writer``, and
 |          ``io.excel.xlsm.writer``.
 |      merge_cells : boolean, default True
 |          Write MultiIndex and Hierarchical Rows as merged cells.
 |      encoding: string, default None
 |          encoding of the resulting excel file. Only necessary for xlwt,
 |          other writers support unicode natively.
 |      inf_rep : string, default 'inf'
 |          Representation for infinity (there is no native representation for
 |          infinity in Excel)
 |      freeze_panes : tuple of integer (length 2), default None
 |          Specifies the one-based bottommost row and rightmost column that
 |          is to be frozen
 |      
 |          .. versionadded:: 0.20.0
 |      
 |      Notes
 |      -----
 |      If passing an existing ExcelWriter object, then the sheet will be added
 |      to the existing workbook.  This can be used to save different
 |      DataFrames to one workbook:
 |      
 |      >>> writer = pd.ExcelWriter('output.xlsx')
 |      >>> df1.to_excel(writer,'Sheet1')
 |      >>> df2.to_excel(writer,'Sheet2')
 |      >>> writer.save()
 |      
 |      For compatibility with to_csv, to_excel serializes lists and dicts to
 |      strings before writing.
 |  
 |  to_feather(self, fname)
 |      write out the binary feather-format for DataFrames
 |      
 |      .. versionadded:: 0.20.0
 |      
 |      Parameters
 |      ----------
 |      fname : str
 |          string file path
 |  
 |  to_gbq(self, destination_table, project_id, chunksize=None, verbose=None, reauth=False, if_exists='fail', private_key=None, auth_local_webserver=False, table_schema=None)
 |      Write a DataFrame to a Google BigQuery table.
 |      
 |      This function requires the `pandas-gbq package
 |      <https://pandas-gbq.readthedocs.io>`__.
 |      
 |      Authentication to the Google BigQuery service is via OAuth 2.0.
 |      
 |      - If ``private_key`` is provided, the library loads the JSON service
 |        account credentials and uses those to authenticate.
 |      
 |      - If no ``private_key`` is provided, the library tries `application
 |        default credentials`_.
 |      
 |        .. _application default credentials:
 |            https://cloud.google.com/docs/authentication/production#providing_credentials_to_your_application
 |      
 |      - If application default credentials are not found or cannot be used
 |        with BigQuery, the library authenticates with user account
 |        credentials. In this case, you will be asked to grant permissions
 |        for product name 'pandas GBQ'.
 |      
 |      Parameters
 |      ----------
 |      destination_table : str
 |          Name of table to be written, in the form 'dataset.tablename'.
 |      project_id : str
 |          Google BigQuery Account project ID.
 |      chunksize : int, optional
 |          Number of rows to be inserted in each chunk from the dataframe.
 |          Set to ``None`` to load the whole dataframe at once.
 |      reauth : bool, default False
 |          Force Google BigQuery to reauthenticate the user. This is useful
 |          if multiple accounts are used.
 |      if_exists : str, default 'fail'
 |          Behavior when the destination table exists. Value can be one of:
 |      
 |          ``'fail'``
 |              If table exists, do nothing.
 |          ``'replace'``
 |              If table exists, drop it, recreate it, and insert data.
 |          ``'append'``
 |              If table exists, insert data. Create if does not exist.
 |      private_key : str, optional
 |          Service account private key in JSON format. Can be file path
 |          or string contents. This is useful for remote server
 |          authentication (eg. Jupyter/IPython notebook on remote host).
 |      auth_local_webserver : bool, default False
 |          Use the `local webserver flow`_ instead of the `console flow`_
 |          when getting user credentials.
 |      
 |          .. _local webserver flow:
 |              http://google-auth-oauthlib.readthedocs.io/en/latest/reference/google_auth_oauthlib.flow.html#google_auth_oauthlib.flow.InstalledAppFlow.run_local_server
 |          .. _console flow:
 |              http://google-auth-oauthlib.readthedocs.io/en/latest/reference/google_auth_oauthlib.flow.html#google_auth_oauthlib.flow.InstalledAppFlow.run_console
 |      
 |          *New in version 0.2.0 of pandas-gbq*.
 |      table_schema : list of dicts, optional
 |          List of BigQuery table fields to which according DataFrame
 |          columns conform to, e.g. ``[{'name': 'col1', 'type':
 |          'STRING'},...]``. If schema is not provided, it will be
 |          generated according to dtypes of DataFrame columns. See
 |          BigQuery API documentation on available names of a field.
 |      
 |          *New in version 0.3.1 of pandas-gbq*.
 |      verbose : boolean, deprecated
 |          *Deprecated in Pandas-GBQ 0.4.0.* Use the `logging module
 |          to adjust verbosity instead
 |          <https://pandas-gbq.readthedocs.io/en/latest/intro.html#logging>`__.
 |      
 |      See Also
 |      --------
 |      pandas_gbq.to_gbq : This function in the pandas-gbq library.
 |      pandas.read_gbq : Read a DataFrame from Google BigQuery.
 |  
 |  to_html(self, buf=None, columns=None, col_space=None, header=True, index=True, na_rep='NaN', formatters=None, float_format=None, sparsify=None, index_names=True, justify=None, bold_rows=True, classes=None, escape=True, max_rows=None, max_cols=None, show_dimensions=False, notebook=False, decimal='.', border=None, table_id=None)
 |      Render a DataFrame as an HTML table.
 |      
 |      `to_html`-specific options:
 |      
 |      bold_rows : boolean, default True
 |          Make the row labels bold in the output
 |      classes : str or list or tuple, default None
 |          CSS class(es) to apply to the resulting html table
 |      escape : boolean, default True
 |          Convert the characters <, >, and & to HTML-safe sequences.
 |      max_rows : int, optional
 |          Maximum number of rows to show before truncating. If None, show
 |          all.
 |      max_cols : int, optional
 |          Maximum number of columns to show before truncating. If None, show
 |          all.
 |      decimal : string, default '.'
 |          Character recognized as decimal separator, e.g. ',' in Europe
 |      
 |          .. versionadded:: 0.18.0
 |      
 |      border : int
 |          A ``border=border`` attribute is included in the opening
 |          `<table>` tag. Default ``pd.options.html.border``.
 |      
 |          .. versionadded:: 0.19.0
 |      
 |      table_id : str, optional
 |          A css id is included in the opening `<table>` tag if specified.
 |      
 |          .. versionadded:: 0.23.0
 |      
 |      
 |      Parameters
 |      ----------
 |      buf : StringIO-like, optional
 |          buffer to write to
 |      columns : sequence, optional
 |          the subset of columns to write; default None writes all columns
 |      col_space : int, optional
 |          the minimum width of each column
 |      header : bool, optional
 |          whether to print column labels, default True
 |      index : bool, optional
 |          whether to print index (row) labels, default True
 |      na_rep : string, optional
 |          string representation of NAN to use, default 'NaN'
 |      formatters : list or dict of one-parameter functions, optional
 |          formatter functions to apply to columns' elements by position or name,
 |          default None. The result of each function must be a unicode string.
 |          List must be of length equal to the number of columns.
 |      float_format : one-parameter function, optional
 |          formatter function to apply to columns' elements if they are floats,
 |          default None. The result of this function must be a unicode string.
 |      sparsify : bool, optional
 |          Set to False for a DataFrame with a hierarchical index to print every
 |          multiindex key at each row, default True
 |      index_names : bool, optional
 |          Prints the names of the indexes, default True
 |      line_width : int, optional
 |          Width to wrap a line in characters, default no wrap
 |      table_id : str, optional
 |          id for the <table> element create by to_html
 |      
 |          .. versionadded:: 0.23.0
 |      justify : str, default None
 |          How to justify the column labels. If None uses the option from
 |          the print configuration (controlled by set_option), 'right' out
 |          of the box. Valid values are
 |      
 |          * left
 |          * right
 |          * center
 |          * justify
 |          * justify-all
 |          * start
 |          * end
 |          * inherit
 |          * match-parent
 |          * initial
 |          * unset
 |      
 |      
 |      Returns
 |      -------
 |      formatted : string (or unicode, depending on data and options)
 |  
 |  to_panel(self)
 |      Transform long (stacked) format (DataFrame) into wide (3D, Panel)
 |      format.
 |      
 |      .. deprecated:: 0.20.0
 |      
 |      Currently the index of the DataFrame must be a 2-level MultiIndex. This
 |      may be generalized later
 |      
 |      Returns
 |      -------
 |      panel : Panel
 |  
 |  to_parquet(self, fname, engine='auto', compression='snappy', **kwargs)
 |      Write a DataFrame to the binary parquet format.
 |      
 |      .. versionadded:: 0.21.0
 |      
 |      This function writes the dataframe as a `parquet file
 |      <https://parquet.apache.org/>`_. You can choose different parquet
 |      backends, and have the option of compression. See
 |      :ref:`the user guide <io.parquet>` for more details.
 |      
 |      Parameters
 |      ----------
 |      fname : str
 |          String file path.
 |      engine : {'auto', 'pyarrow', 'fastparquet'}, default 'auto'
 |          Parquet library to use. If 'auto', then the option
 |          ``io.parquet.engine`` is used. The default ``io.parquet.engine``
 |          behavior is to try 'pyarrow', falling back to 'fastparquet' if
 |          'pyarrow' is unavailable.
 |      compression : {'snappy', 'gzip', 'brotli', None}, default 'snappy'
 |          Name of the compression to use. Use ``None`` for no compression.
 |      **kwargs
 |          Additional arguments passed to the parquet library. See
 |          :ref:`pandas io <io.parquet>` for more details.
 |      
 |      See Also
 |      --------
 |      read_parquet : Read a parquet file.
 |      DataFrame.to_csv : Write a csv file.
 |      DataFrame.to_sql : Write to a sql table.
 |      DataFrame.to_hdf : Write to hdf.
 |      
 |      Notes
 |      -----
 |      This function requires either the `fastparquet
 |      <https://pypi.org/project/fastparquet>`_ or `pyarrow
 |      <https://arrow.apache.org/docs/python/>`_ library.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame(data={'col1': [1, 2], 'col2': [3, 4]})
 |      >>> df.to_parquet('df.parquet.gzip', compression='gzip')
 |      >>> pd.read_parquet('df.parquet.gzip')
 |         col1  col2
 |      0     1     3
 |      1     2     4
 |  
 |  to_period(self, freq=None, axis=0, copy=True)
 |      Convert DataFrame from DatetimeIndex to PeriodIndex with desired
 |      frequency (inferred from index if not passed)
 |      
 |      Parameters
 |      ----------
 |      freq : string, default
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          The axis to convert (the index by default)
 |      copy : boolean, default True
 |          If False then underlying input data is not copied
 |      
 |      Returns
 |      -------
 |      ts : TimeSeries with PeriodIndex
 |  
 |  to_records(self, index=True, convert_datetime64=None)
 |      Convert DataFrame to a NumPy record array.
 |      
 |      Index will be put in the 'index' field of the record array if
 |      requested.
 |      
 |      Parameters
 |      ----------
 |      index : boolean, default True
 |          Include index in resulting record array, stored in 'index' field.
 |      convert_datetime64 : boolean, default None
 |          .. deprecated:: 0.23.0
 |      
 |          Whether to convert the index to datetime.datetime if it is a
 |          DatetimeIndex.
 |      
 |      Returns
 |      -------
 |      y : numpy.recarray
 |      
 |      See Also
 |      --------
 |      DataFrame.from_records: convert structured or record ndarray
 |          to DataFrame.
 |      numpy.recarray: ndarray that allows field access using
 |          attributes, analogous to typed columns in a
 |          spreadsheet.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'A': [1, 2], 'B': [0.5, 0.75]},
 |      ...                   index=['a', 'b'])
 |      >>> df
 |         A     B
 |      a  1  0.50
 |      b  2  0.75
 |      >>> df.to_records()
 |      rec.array([('a', 1, 0.5 ), ('b', 2, 0.75)],
 |                dtype=[('index', 'O'), ('A', '<i8'), ('B', '<f8')])
 |      
 |      The index can be excluded from the record array:
 |      
 |      >>> df.to_records(index=False)
 |      rec.array([(1, 0.5 ), (2, 0.75)],
 |                dtype=[('A', '<i8'), ('B', '<f8')])
 |      
 |      By default, timestamps are converted to `datetime.datetime`:
 |      
 |      >>> df.index = pd.date_range('2018-01-01 09:00', periods=2, freq='min')
 |      >>> df
 |                           A     B
 |      2018-01-01 09:00:00  1  0.50
 |      2018-01-01 09:01:00  2  0.75
 |      >>> df.to_records()
 |      rec.array([(datetime.datetime(2018, 1, 1, 9, 0), 1, 0.5 ),
 |                 (datetime.datetime(2018, 1, 1, 9, 1), 2, 0.75)],
 |                dtype=[('index', 'O'), ('A', '<i8'), ('B', '<f8')])
 |      
 |      The timestamp conversion can be disabled so NumPy's datetime64
 |      data type is used instead:
 |      
 |      >>> df.to_records(convert_datetime64=False)
 |      rec.array([('2018-01-01T09:00:00.000000000', 1, 0.5 ),
 |                 ('2018-01-01T09:01:00.000000000', 2, 0.75)],
 |                dtype=[('index', '<M8[ns]'), ('A', '<i8'), ('B', '<f8')])
 |  
 |  to_sparse(self, fill_value=None, kind='block')
 |      Convert to SparseDataFrame
 |      
 |      Parameters
 |      ----------
 |      fill_value : float, default NaN
 |      kind : {'block', 'integer'}
 |      
 |      Returns
 |      -------
 |      y : SparseDataFrame
 |  
 |  to_stata(self, fname, convert_dates=None, write_index=True, encoding='latin-1', byteorder=None, time_stamp=None, data_label=None, variable_labels=None, version=114, convert_strl=None)
 |      Export Stata binary dta files.
 |      
 |      Parameters
 |      ----------
 |      fname : path (string), buffer or path object
 |          string, path object (pathlib.Path or py._path.local.LocalPath) or
 |          object implementing a binary write() functions. If using a buffer
 |          then the buffer will not be automatically closed after the file
 |          data has been written.
 |      convert_dates : dict
 |          Dictionary mapping columns containing datetime types to stata
 |          internal format to use when writing the dates. Options are 'tc',
 |          'td', 'tm', 'tw', 'th', 'tq', 'ty'. Column can be either an integer
 |          or a name. Datetime columns that do not have a conversion type
 |          specified will be converted to 'tc'. Raises NotImplementedError if
 |          a datetime column has timezone information.
 |      write_index : bool
 |          Write the index to Stata dataset.
 |      encoding : str
 |          Default is latin-1. Unicode is not supported.
 |      byteorder : str
 |          Can be ">", "<", "little", or "big". default is `sys.byteorder`.
 |      time_stamp : datetime
 |          A datetime to use as file creation date.  Default is the current
 |          time.
 |      data_label : str
 |          A label for the data set.  Must be 80 characters or smaller.
 |      variable_labels : dict
 |          Dictionary containing columns as keys and variable labels as
 |          values. Each label must be 80 characters or smaller.
 |      
 |          .. versionadded:: 0.19.0
 |      
 |      version : {114, 117}
 |          Version to use in the output dta file.  Version 114 can be used
 |          read by Stata 10 and later.  Version 117 can be read by Stata 13
 |          or later. Version 114 limits string variables to 244 characters or
 |          fewer while 117 allows strings with lengths up to 2,000,000
 |          characters.
 |      
 |          .. versionadded:: 0.23.0
 |      
 |      convert_strl : list, optional
 |          List of column names to convert to string columns to Stata StrL
 |          format. Only available if version is 117.  Storing strings in the
 |          StrL format can produce smaller dta files if strings have more than
 |          8 characters and values are repeated.
 |      
 |          .. versionadded:: 0.23.0
 |      
 |      Raises
 |      ------
 |      NotImplementedError
 |          * If datetimes contain timezone information
 |          * Column dtype is not representable in Stata
 |      ValueError
 |          * Columns listed in convert_dates are neither datetime64[ns]
 |            or datetime.datetime
 |          * Column listed in convert_dates is not in DataFrame
 |          * Categorical label contains more than 32,000 characters
 |      
 |          .. versionadded:: 0.19.0
 |      
 |      See Also
 |      --------
 |      pandas.read_stata : Import Stata data files
 |      pandas.io.stata.StataWriter : low-level writer for Stata data files
 |      pandas.io.stata.StataWriter117 : low-level writer for version 117 files
 |      
 |      Examples
 |      --------
 |      >>> data.to_stata('./data_file.dta')
 |      
 |      Or with dates
 |      
 |      >>> data.to_stata('./date_data_file.dta', {2 : 'tw'})
 |      
 |      Alternatively you can create an instance of the StataWriter class
 |      
 |      >>> writer = StataWriter('./data_file.dta', data)
 |      >>> writer.write_file()
 |      
 |      With dates:
 |      
 |      >>> writer = StataWriter('./date_data_file.dta', data, {2 : 'tw'})
 |      >>> writer.write_file()
 |  
 |  to_string(self, buf=None, columns=None, col_space=None, header=True, index=True, na_rep='NaN', formatters=None, float_format=None, sparsify=None, index_names=True, justify=None, line_width=None, max_rows=None, max_cols=None, show_dimensions=False)
 |      Render a DataFrame to a console-friendly tabular output.
 |      
 |      Parameters
 |      ----------
 |      buf : StringIO-like, optional
 |          buffer to write to
 |      columns : sequence, optional
 |          the subset of columns to write; default None writes all columns
 |      col_space : int, optional
 |          the minimum width of each column
 |      header : bool, optional
 |          Write out the column names. If a list of strings is given, it is assumed to be aliases for the column names
 |      index : bool, optional
 |          whether to print index (row) labels, default True
 |      na_rep : string, optional
 |          string representation of NAN to use, default 'NaN'
 |      formatters : list or dict of one-parameter functions, optional
 |          formatter functions to apply to columns' elements by position or name,
 |          default None. The result of each function must be a unicode string.
 |          List must be of length equal to the number of columns.
 |      float_format : one-parameter function, optional
 |          formatter function to apply to columns' elements if they are floats,
 |          default None. The result of this function must be a unicode string.
 |      sparsify : bool, optional
 |          Set to False for a DataFrame with a hierarchical index to print every
 |          multiindex key at each row, default True
 |      index_names : bool, optional
 |          Prints the names of the indexes, default True
 |      line_width : int, optional
 |          Width to wrap a line in characters, default no wrap
 |      table_id : str, optional
 |          id for the <table> element create by to_html
 |      
 |          .. versionadded:: 0.23.0
 |      justify : str, default None
 |          How to justify the column labels. If None uses the option from
 |          the print configuration (controlled by set_option), 'right' out
 |          of the box. Valid values are
 |      
 |          * left
 |          * right
 |          * center
 |          * justify
 |          * justify-all
 |          * start
 |          * end
 |          * inherit
 |          * match-parent
 |          * initial
 |          * unset
 |      
 |      
 |      Returns
 |      -------
 |      formatted : string (or unicode, depending on data and options)
 |  
 |  to_timestamp(self, freq=None, how='start', axis=0, copy=True)
 |      Cast to DatetimeIndex of timestamps, at *beginning* of period
 |      
 |      Parameters
 |      ----------
 |      freq : string, default frequency of PeriodIndex
 |          Desired frequency
 |      how : {'s', 'e', 'start', 'end'}
 |          Convention for converting period to timestamp; start of period
 |          vs. end
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          The axis to convert (the index by default)
 |      copy : boolean, default True
 |          If false then underlying input data is not copied
 |      
 |      Returns
 |      -------
 |      df : DataFrame with DatetimeIndex
 |  
 |  transform(self, func, *args, **kwargs)
 |      Call function producing a like-indexed NDFrame
 |      and return a NDFrame with the transformed values
 |      
 |      .. versionadded:: 0.20.0
 |      
 |      Parameters
 |      ----------
 |      func : callable, string, dictionary, or list of string/callables
 |          To apply to column
 |      
 |          Accepted Combinations are:
 |      
 |          - string function name
 |          - function
 |          - list of functions
 |          - dict of column names -> functions (or list of functions)
 |      
 |      Returns
 |      -------
 |      transformed : NDFrame
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame(np.random.randn(10, 3), columns=['A', 'B', 'C'],
 |      ...                   index=pd.date_range('1/1/2000', periods=10))
 |      df.iloc[3:7] = np.nan
 |      
 |      >>> df.transform(lambda x: (x - x.mean()) / x.std())
 |                         A         B         C
 |      2000-01-01  0.579457  1.236184  0.123424
 |      2000-01-02  0.370357 -0.605875 -1.231325
 |      2000-01-03  1.455756 -0.277446  0.288967
 |      2000-01-04       NaN       NaN       NaN
 |      2000-01-05       NaN       NaN       NaN
 |      2000-01-06       NaN       NaN       NaN
 |      2000-01-07       NaN       NaN       NaN
 |      2000-01-08 -0.498658  1.274522  1.642524
 |      2000-01-09 -0.540524 -1.012676 -0.828968
 |      2000-01-10 -1.366388 -0.614710  0.005378
 |      
 |      See also
 |      --------
 |      pandas.NDFrame.aggregate
 |      pandas.NDFrame.apply
 |  
 |  transpose(self, *args, **kwargs)
 |      Transpose index and columns.
 |      
 |      Reflect the DataFrame over its main diagonal by writing rows as columns
 |      and vice-versa. The property :attr:`.T` is an accessor to the method
 |      :meth:`transpose`.
 |      
 |      Parameters
 |      ----------
 |      copy : bool, default False
 |          If True, the underlying data is copied. Otherwise (default), no
 |          copy is made if possible.
 |      *args, **kwargs
 |          Additional keywords have no effect but might be accepted for
 |          compatibility with numpy.
 |      
 |      Returns
 |      -------
 |      DataFrame
 |          The transposed DataFrame.
 |      
 |      See Also
 |      --------
 |      numpy.transpose : Permute the dimensions of a given array.
 |      
 |      Notes
 |      -----
 |      Transposing a DataFrame with mixed dtypes will result in a homogeneous
 |      DataFrame with the `object` dtype. In such a case, a copy of the data
 |      is always made.
 |      
 |      Examples
 |      --------
 |      **Square DataFrame with homogeneous dtype**
 |      
 |      >>> d1 = {'col1': [1, 2], 'col2': [3, 4]}
 |      >>> df1 = pd.DataFrame(data=d1)
 |      >>> df1
 |         col1  col2
 |      0     1     3
 |      1     2     4
 |      
 |      >>> df1_transposed = df1.T # or df1.transpose()
 |      >>> df1_transposed
 |            0  1
 |      col1  1  2
 |      col2  3  4
 |      
 |      When the dtype is homogeneous in the original DataFrame, we get a
 |      transposed DataFrame with the same dtype:
 |      
 |      >>> df1.dtypes
 |      col1    int64
 |      col2    int64
 |      dtype: object
 |      >>> df1_transposed.dtypes
 |      0    int64
 |      1    int64
 |      dtype: object
 |      
 |      **Non-square DataFrame with mixed dtypes**
 |      
 |      >>> d2 = {'name': ['Alice', 'Bob'],
 |      ...       'score': [9.5, 8],
 |      ...       'employed': [False, True],
 |      ...       'kids': [0, 0]}
 |      >>> df2 = pd.DataFrame(data=d2)
 |      >>> df2
 |          name  score  employed  kids
 |      0  Alice    9.5     False     0
 |      1    Bob    8.0      True     0
 |      
 |      >>> df2_transposed = df2.T # or df2.transpose()
 |      >>> df2_transposed
 |                    0     1
 |      name      Alice   Bob
 |      score       9.5     8
 |      employed  False  True
 |      kids          0     0
 |      
 |      When the DataFrame has mixed dtypes, we get a transposed DataFrame with
 |      the `object` dtype:
 |      
 |      >>> df2.dtypes
 |      name         object
 |      score       float64
 |      employed       bool
 |      kids          int64
 |      dtype: object
 |      >>> df2_transposed.dtypes
 |      0    object
 |      1    object
 |      dtype: object
 |  
 |  truediv(self, other, axis='columns', level=None, fill_value=None)
 |      Floating division of dataframe and other, element-wise (binary operator `truediv`).
 |      
 |      Equivalent to ``dataframe / other``, but with support to substitute a fill_value for
 |      missing data in one of the inputs.
 |      
 |      Parameters
 |      ----------
 |      other : Series, DataFrame, or constant
 |      axis : {0, 1, 'index', 'columns'}
 |          For Series input, axis to match Series index on
 |      level : int or name
 |          Broadcast across a level, matching Index values on the
 |          passed MultiIndex level
 |      fill_value : None or float value, default None
 |          Fill existing missing (NaN) values, and any new element needed for
 |          successful DataFrame alignment, with this value before computation.
 |          If data in both corresponding DataFrame locations is missing
 |          the result will be missing
 |      
 |      Notes
 |      -----
 |      Mismatched indices will be unioned together
 |      
 |      Returns
 |      -------
 |      result : DataFrame
 |      
 |      Examples
 |      --------
 |      None
 |      
 |      See also
 |      --------
 |      DataFrame.rtruediv
 |  
 |  unstack(self, level=-1, fill_value=None)
 |      Pivot a level of the (necessarily hierarchical) index labels, returning
 |      a DataFrame having a new level of column labels whose inner-most level
 |      consists of the pivoted index labels. If the index is not a MultiIndex,
 |      the output will be a Series (the analogue of stack when the columns are
 |      not a MultiIndex).
 |      The level involved will automatically get sorted.
 |      
 |      Parameters
 |      ----------
 |      level : int, string, or list of these, default -1 (last level)
 |          Level(s) of index to unstack, can pass level name
 |      fill_value : replace NaN with this value if the unstack produces
 |          missing values
 |      
 |          .. versionadded:: 0.18.0
 |      
 |      See also
 |      --------
 |      DataFrame.pivot : Pivot a table based on column values.
 |      DataFrame.stack : Pivot a level of the column labels (inverse operation
 |          from `unstack`).
 |      
 |      Examples
 |      --------
 |      >>> index = pd.MultiIndex.from_tuples([('one', 'a'), ('one', 'b'),
 |      ...                                    ('two', 'a'), ('two', 'b')])
 |      >>> s = pd.Series(np.arange(1.0, 5.0), index=index)
 |      >>> s
 |      one  a   1.0
 |           b   2.0
 |      two  a   3.0
 |           b   4.0
 |      dtype: float64
 |      
 |      >>> s.unstack(level=-1)
 |           a   b
 |      one  1.0  2.0
 |      two  3.0  4.0
 |      
 |      >>> s.unstack(level=0)
 |         one  two
 |      a  1.0   3.0
 |      b  2.0   4.0
 |      
 |      >>> df = s.unstack(level=0)
 |      >>> df.unstack()
 |      one  a  1.0
 |           b  2.0
 |      two  a  3.0
 |           b  4.0
 |      dtype: float64
 |      
 |      Returns
 |      -------
 |      unstacked : DataFrame or Series
 |  
 |  update(self, other, join='left', overwrite=True, filter_func=None, raise_conflict=False)
 |      Modify in place using non-NA values from another DataFrame.
 |      
 |      Aligns on indices. There is no return value.
 |      
 |      Parameters
 |      ----------
 |      other : DataFrame, or object coercible into a DataFrame
 |          Should have at least one matching index/column label
 |          with the original DataFrame. If a Series is passed,
 |          its name attribute must be set, and that will be
 |          used as the column name to align with the original DataFrame.
 |      join : {'left'}, default 'left'
 |          Only left join is implemented, keeping the index and columns of the
 |          original object.
 |      overwrite : bool, default True
 |          How to handle non-NA values for overlapping keys:
 |      
 |          * True: overwrite original DataFrame's values
 |            with values from `other`.
 |          * False: only update values that are NA in
 |            the original DataFrame.
 |      
 |      filter_func : callable(1d-array) -> boolean 1d-array, optional
 |          Can choose to replace values other than NA. Return True for values
 |          that should be updated.
 |      raise_conflict : bool, default False
 |          If True, will raise a ValueError if the DataFrame and `other`
 |          both contain non-NA data in the same place.
 |      
 |      Raises
 |      ------
 |      ValueError
 |          When `raise_conflict` is True and there's overlapping non-NA data.
 |      
 |      See Also
 |      --------
 |      dict.update : Similar method for dictionaries.
 |      DataFrame.merge : For column(s)-on-columns(s) operations.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'A': [1, 2, 3],
 |      ...                    'B': [400, 500, 600]})
 |      >>> new_df = pd.DataFrame({'B': [4, 5, 6],
 |      ...                        'C': [7, 8, 9]})
 |      >>> df.update(new_df)
 |      >>> df
 |         A  B
 |      0  1  4
 |      1  2  5
 |      2  3  6
 |      
 |      The DataFrame's length does not increase as a result of the update,
 |      only values at matching index/column labels are updated.
 |      
 |      >>> df = pd.DataFrame({'A': ['a', 'b', 'c'],
 |      ...                    'B': ['x', 'y', 'z']})
 |      >>> new_df = pd.DataFrame({'B': ['d', 'e', 'f', 'g', 'h', 'i']})
 |      >>> df.update(new_df)
 |      >>> df
 |         A  B
 |      0  a  d
 |      1  b  e
 |      2  c  f
 |      
 |      For Series, it's name attribute must be set.
 |      
 |      >>> df = pd.DataFrame({'A': ['a', 'b', 'c'],
 |      ...                    'B': ['x', 'y', 'z']})
 |      >>> new_column = pd.Series(['d', 'e'], name='B', index=[0, 2])
 |      >>> df.update(new_column)
 |      >>> df
 |         A  B
 |      0  a  d
 |      1  b  y
 |      2  c  e
 |      >>> df = pd.DataFrame({'A': ['a', 'b', 'c'],
 |      ...                    'B': ['x', 'y', 'z']})
 |      >>> new_df = pd.DataFrame({'B': ['d', 'e']}, index=[1, 2])
 |      >>> df.update(new_df)
 |      >>> df
 |         A  B
 |      0  a  x
 |      1  b  d
 |      2  c  e
 |      
 |      If `other` contains NaNs the corresponding values are not updated
 |      in the original dataframe.
 |      
 |      >>> df = pd.DataFrame({'A': [1, 2, 3],
 |      ...                    'B': [400, 500, 600]})
 |      >>> new_df = pd.DataFrame({'B': [4, np.nan, 6]})
 |      >>> df.update(new_df)
 |      >>> df
 |         A      B
 |      0  1    4.0
 |      1  2  500.0
 |      2  3    6.0
 |  
 |  var(self, axis=None, skipna=None, level=None, ddof=1, numeric_only=None, **kwargs)
 |      Return unbiased variance over requested axis.
 |      
 |      Normalized by N-1 by default. This can be changed using the ddof argument
 |      
 |      Parameters
 |      ----------
 |      axis : {index (0), columns (1)}
 |      skipna : boolean, default True
 |          Exclude NA/null values. If an entire row/column is NA, the result
 |          will be NA
 |      level : int or level name, default None
 |          If the axis is a MultiIndex (hierarchical), count along a
 |          particular level, collapsing into a Series
 |      ddof : int, default 1
 |          Delta Degrees of Freedom. The divisor used in calculations is N - ddof,
 |          where N represents the number of elements.
 |      numeric_only : boolean, default None
 |          Include only float, int, boolean columns. If None, will attempt to use
 |          everything, then use only numeric data. Not implemented for Series.
 |      
 |      Returns
 |      -------
 |      var : Series or DataFrame (if level specified)
 |  
 |  ----------------------------------------------------------------------
 |  Class methods defined here:
 |  
 |  from_csv(path, header=0, sep=',', index_col=0, parse_dates=True, encoding=None, tupleize_cols=None, infer_datetime_format=False) from builtins.type
 |      Read CSV file.
 |      
 |      .. deprecated:: 0.21.0
 |          Use :func:`pandas.read_csv` instead.
 |      
 |      It is preferable to use the more powerful :func:`pandas.read_csv`
 |      for most general purposes, but ``from_csv`` makes for an easy
 |      roundtrip to and from a file (the exact counterpart of
 |      ``to_csv``), especially with a DataFrame of time series data.
 |      
 |      This method only differs from the preferred :func:`pandas.read_csv`
 |      in some defaults:
 |      
 |      - `index_col` is ``0`` instead of ``None`` (take first column as index
 |        by default)
 |      - `parse_dates` is ``True`` instead of ``False`` (try parsing the index
 |        as datetime by default)
 |      
 |      So a ``pd.DataFrame.from_csv(path)`` can be replaced by
 |      ``pd.read_csv(path, index_col=0, parse_dates=True)``.
 |      
 |      Parameters
 |      ----------
 |      path : string file path or file handle / StringIO
 |      header : int, default 0
 |          Row to use as header (skip prior rows)
 |      sep : string, default ','
 |          Field delimiter
 |      index_col : int or sequence, default 0
 |          Column to use for index. If a sequence is given, a MultiIndex
 |          is used. Different default from read_table
 |      parse_dates : boolean, default True
 |          Parse dates. Different default from read_table
 |      tupleize_cols : boolean, default False
 |          write multi_index columns as a list of tuples (if True)
 |          or new (expanded format) if False)
 |      infer_datetime_format: boolean, default False
 |          If True and `parse_dates` is True for a column, try to infer the
 |          datetime format based on the first datetime string. If the format
 |          can be inferred, there often will be a large parsing speed-up.
 |      
 |      See also
 |      --------
 |      pandas.read_csv
 |      
 |      Returns
 |      -------
 |      y : DataFrame
 |  
 |  from_dict(data, orient='columns', dtype=None, columns=None) from builtins.type
 |      Construct DataFrame from dict of array-like or dicts.
 |      
 |      Creates DataFrame object from dictionary by columns or by index
 |      allowing dtype specification.
 |      
 |      Parameters
 |      ----------
 |      data : dict
 |          Of the form {field : array-like} or {field : dict}.
 |      orient : {'columns', 'index'}, default 'columns'
 |          The "orientation" of the data. If the keys of the passed dict
 |          should be the columns of the resulting DataFrame, pass 'columns'
 |          (default). Otherwise if the keys should be rows, pass 'index'.
 |      dtype : dtype, default None
 |          Data type to force, otherwise infer.
 |      columns : list, default None
 |          Column labels to use when ``orient='index'``. Raises a ValueError
 |          if used with ``orient='columns'``.
 |      
 |          .. versionadded:: 0.23.0
 |      
 |      Returns
 |      -------
 |      pandas.DataFrame
 |      
 |      See Also
 |      --------
 |      DataFrame.from_records : DataFrame from ndarray (structured
 |          dtype), list of tuples, dict, or DataFrame
 |      DataFrame : DataFrame object creation using constructor
 |      
 |      Examples
 |      --------
 |      By default the keys of the dict become the DataFrame columns:
 |      
 |      >>> data = {'col_1': [3, 2, 1, 0], 'col_2': ['a', 'b', 'c', 'd']}
 |      >>> pd.DataFrame.from_dict(data)
 |         col_1 col_2
 |      0      3     a
 |      1      2     b
 |      2      1     c
 |      3      0     d
 |      
 |      Specify ``orient='index'`` to create the DataFrame using dictionary
 |      keys as rows:
 |      
 |      >>> data = {'row_1': [3, 2, 1, 0], 'row_2': ['a', 'b', 'c', 'd']}
 |      >>> pd.DataFrame.from_dict(data, orient='index')
 |             0  1  2  3
 |      row_1  3  2  1  0
 |      row_2  a  b  c  d
 |      
 |      When using the 'index' orientation, the column names can be
 |      specified manually:
 |      
 |      >>> pd.DataFrame.from_dict(data, orient='index',
 |      ...                        columns=['A', 'B', 'C', 'D'])
 |             A  B  C  D
 |      row_1  3  2  1  0
 |      row_2  a  b  c  d
 |  
 |  from_items(items, columns=None, orient='columns') from builtins.type
 |      Construct a dataframe from a list of tuples
 |      
 |      .. deprecated:: 0.23.0
 |        `from_items` is deprecated and will be removed in a future version.
 |        Use :meth:`DataFrame.from_dict(dict(items)) <DataFrame.from_dict>`
 |        instead.
 |        :meth:`DataFrame.from_dict(OrderedDict(items)) <DataFrame.from_dict>`
 |        may be used to preserve the key order.
 |      
 |      Convert (key, value) pairs to DataFrame. The keys will be the axis
 |      index (usually the columns, but depends on the specified
 |      orientation). The values should be arrays or Series.
 |      
 |      Parameters
 |      ----------
 |      items : sequence of (key, value) pairs
 |          Values should be arrays or Series.
 |      columns : sequence of column labels, optional
 |          Must be passed if orient='index'.
 |      orient : {'columns', 'index'}, default 'columns'
 |          The "orientation" of the data. If the keys of the
 |          input correspond to column labels, pass 'columns'
 |          (default). Otherwise if the keys correspond to the index,
 |          pass 'index'.
 |      
 |      Returns
 |      -------
 |      frame : DataFrame
 |  
 |  from_records(data, index=None, exclude=None, columns=None, coerce_float=False, nrows=None) from builtins.type
 |      Convert structured or record ndarray to DataFrame
 |      
 |      Parameters
 |      ----------
 |      data : ndarray (structured dtype), list of tuples, dict, or DataFrame
 |      index : string, list of fields, array-like
 |          Field of array to use as the index, alternately a specific set of
 |          input labels to use
 |      exclude : sequence, default None
 |          Columns or fields to exclude
 |      columns : sequence, default None
 |          Column names to use. If the passed data do not have names
 |          associated with them, this argument provides names for the
 |          columns. Otherwise this argument indicates the order of the columns
 |          in the result (any names not found in the data will become all-NA
 |          columns)
 |      coerce_float : boolean, default False
 |          Attempt to convert values of non-string, non-numeric objects (like
 |          decimal.Decimal) to floating point, useful for SQL result sets
 |      
 |      Returns
 |      -------
 |      df : DataFrame
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  T
 |      Transpose index and columns.
 |      
 |      Reflect the DataFrame over its main diagonal by writing rows as columns
 |      and vice-versa. The property :attr:`.T` is an accessor to the method
 |      :meth:`transpose`.
 |      
 |      Parameters
 |      ----------
 |      copy : bool, default False
 |          If True, the underlying data is copied. Otherwise (default), no
 |          copy is made if possible.
 |      *args, **kwargs
 |          Additional keywords have no effect but might be accepted for
 |          compatibility with numpy.
 |      
 |      Returns
 |      -------
 |      DataFrame
 |          The transposed DataFrame.
 |      
 |      See Also
 |      --------
 |      numpy.transpose : Permute the dimensions of a given array.
 |      
 |      Notes
 |      -----
 |      Transposing a DataFrame with mixed dtypes will result in a homogeneous
 |      DataFrame with the `object` dtype. In such a case, a copy of the data
 |      is always made.
 |      
 |      Examples
 |      --------
 |      **Square DataFrame with homogeneous dtype**
 |      
 |      >>> d1 = {'col1': [1, 2], 'col2': [3, 4]}
 |      >>> df1 = pd.DataFrame(data=d1)
 |      >>> df1
 |         col1  col2
 |      0     1     3
 |      1     2     4
 |      
 |      >>> df1_transposed = df1.T # or df1.transpose()
 |      >>> df1_transposed
 |            0  1
 |      col1  1  2
 |      col2  3  4
 |      
 |      When the dtype is homogeneous in the original DataFrame, we get a
 |      transposed DataFrame with the same dtype:
 |      
 |      >>> df1.dtypes
 |      col1    int64
 |      col2    int64
 |      dtype: object
 |      >>> df1_transposed.dtypes
 |      0    int64
 |      1    int64
 |      dtype: object
 |      
 |      **Non-square DataFrame with mixed dtypes**
 |      
 |      >>> d2 = {'name': ['Alice', 'Bob'],
 |      ...       'score': [9.5, 8],
 |      ...       'employed': [False, True],
 |      ...       'kids': [0, 0]}
 |      >>> df2 = pd.DataFrame(data=d2)
 |      >>> df2
 |          name  score  employed  kids
 |      0  Alice    9.5     False     0
 |      1    Bob    8.0      True     0
 |      
 |      >>> df2_transposed = df2.T # or df2.transpose()
 |      >>> df2_transposed
 |                    0     1
 |      name      Alice   Bob
 |      score       9.5     8
 |      employed  False  True
 |      kids          0     0
 |      
 |      When the DataFrame has mixed dtypes, we get a transposed DataFrame with
 |      the `object` dtype:
 |      
 |      >>> df2.dtypes
 |      name         object
 |      score       float64
 |      employed       bool
 |      kids          int64
 |      dtype: object
 |      >>> df2_transposed.dtypes
 |      0    object
 |      1    object
 |      dtype: object
 |  
 |  axes
 |      Return a list representing the axes of the DataFrame.
 |      
 |      It has the row axis labels and column axis labels as the only members.
 |      They are returned in that order.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
 |      >>> df.axes
 |      [RangeIndex(start=0, stop=2, step=1), Index(['coll', 'col2'],
 |      dtype='object')]
 |  
 |  columns
 |      The column labels of the DataFrame.
 |  
 |  index
 |      The index (row labels) of the DataFrame.
 |  
 |  shape
 |      Return a tuple representing the dimensionality of the DataFrame.
 |      
 |      See Also
 |      --------
 |      ndarray.shape
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
 |      >>> df.shape
 |      (2, 2)
 |      
 |      >>> df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4],
 |      ...                    'col3': [5, 6]})
 |      >>> df.shape
 |      (2, 3)
 |  
 |  style
 |      Property returning a Styler object containing methods for
 |      building a styled HTML representation fo the DataFrame.
 |      
 |      See Also
 |      --------
 |      pandas.io.formats.style.Styler
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  plot = <class 'pandas.plotting._core.FramePlotMethods'>
 |      DataFrame plotting accessor and method
 |      
 |      Examples
 |      --------
 |      >>> df.plot.line()
 |      >>> df.plot.scatter('x', 'y')
 |      >>> df.plot.hexbin()
 |      
 |      These plotting methods can also be accessed by calling the accessor as a
 |      method with the ``kind`` argument:
 |      ``df.plot(kind='line')`` is equivalent to ``df.plot.line()``
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from pandas.core.generic.NDFrame:
 |  
 |  __abs__(self)
 |  
 |  __array__(self, dtype=None)
 |  
 |  __array_wrap__(self, result, context=None)
 |  
 |  __bool__ = __nonzero__(self)
 |  
 |  __contains__(self, key)
 |      True if the key is in the info axis
 |  
 |  __copy__(self, deep=True)
 |  
 |  __deepcopy__(self, memo=None)
 |  
 |  __delitem__(self, key)
 |      Delete item
 |  
 |  __finalize__(self, other, method=None, **kwargs)
 |      Propagate metadata from other to self.
 |      
 |      Parameters
 |      ----------
 |      other : the object from which to get the attributes that we are going
 |          to propagate
 |      method : optional, a passed method name ; possibly to take different
 |          types of propagation actions based on this
 |  
 |  __getattr__(self, name)
 |      After regular attribute access, try looking up the name
 |      This allows simpler access to columns for interactive use.
 |  
 |  __getstate__(self)
 |  
 |  __hash__(self)
 |      Return hash(self).
 |  
 |  __invert__(self)
 |  
 |  __iter__(self)
 |      Iterate over infor axis
 |  
 |  __neg__(self)
 |  
 |  __nonzero__(self)
 |  
 |  __pos__(self)
 |  
 |  __round__(self, decimals=0)
 |  
 |  __setattr__(self, name, value)
 |      After regular attribute access, try setting the name
 |      This allows simpler access to columns for interactive use.
 |  
 |  __setstate__(self, state)
 |  
 |  abs(self)
 |      Return a Series/DataFrame with absolute numeric value of each element.
 |      
 |      This function only applies to elements that are all numeric.
 |      
 |      Returns
 |      -------
 |      abs
 |          Series/DataFrame containing the absolute value of each element.
 |      
 |      Notes
 |      -----
 |      For ``complex`` inputs, ``1.2 + 1j``, the absolute value is
 |      :math:`\sqrt{ a^2 + b^2 }`.
 |      
 |      Examples
 |      --------
 |      Absolute numeric values in a Series.
 |      
 |      >>> s = pd.Series([-1.10, 2, -3.33, 4])
 |      >>> s.abs()
 |      0    1.10
 |      1    2.00
 |      2    3.33
 |      3    4.00
 |      dtype: float64
 |      
 |      Absolute numeric values in a Series with complex numbers.
 |      
 |      >>> s = pd.Series([1.2 + 1j])
 |      >>> s.abs()
 |      0    1.56205
 |      dtype: float64
 |      
 |      Absolute numeric values in a Series with a Timedelta element.
 |      
 |      >>> s = pd.Series([pd.Timedelta('1 days')])
 |      >>> s.abs()
 |      0   1 days
 |      dtype: timedelta64[ns]
 |      
 |      Select rows with data closest to certain value using argsort (from
 |      `StackOverflow <https://stackoverflow.com/a/17758115>`__).
 |      
 |      >>> df = pd.DataFrame({
 |      ...     'a': [4, 5, 6, 7],
 |      ...     'b': [10, 20, 30, 40],
 |      ...     'c': [100, 50, -30, -50]
 |      ... })
 |      >>> df
 |           a    b    c
 |      0    4   10  100
 |      1    5   20   50
 |      2    6   30  -30
 |      3    7   40  -50
 |      >>> df.loc[(df.c - 43).abs().argsort()]
 |           a    b    c
 |      1    5   20   50
 |      0    4   10  100
 |      2    6   30  -30
 |      3    7   40  -50
 |      
 |      See Also
 |      --------
 |      numpy.absolute : calculate the absolute value element-wise.
 |  
 |  add_prefix(self, prefix)
 |      Prefix labels with string `prefix`.
 |      
 |      For Series, the row labels are prefixed.
 |      For DataFrame, the column labels are prefixed.
 |      
 |      Parameters
 |      ----------
 |      prefix : str
 |          The string to add before each label.
 |      
 |      Returns
 |      -------
 |      Series or DataFrame
 |          New Series or DataFrame with updated labels.
 |      
 |      See Also
 |      --------
 |      Series.add_suffix: Suffix row labels with string `suffix`.
 |      DataFrame.add_suffix: Suffix column labels with string `suffix`.
 |      
 |      Examples
 |      --------
 |      >>> s = pd.Series([1, 2, 3, 4])
 |      >>> s
 |      0    1
 |      1    2
 |      2    3
 |      3    4
 |      dtype: int64
 |      
 |      >>> s.add_prefix('item_')
 |      item_0    1
 |      item_1    2
 |      item_2    3
 |      item_3    4
 |      dtype: int64
 |      
 |      >>> df = pd.DataFrame({'A': [1, 2, 3, 4],  'B': [3, 4, 5, 6]})
 |      >>> df
 |         A  B
 |      0  1  3
 |      1  2  4
 |      2  3  5
 |      3  4  6
 |      
 |      >>> df.add_prefix('col_')
 |           col_A  col_B
 |      0       1       3
 |      1       2       4
 |      2       3       5
 |      3       4       6
 |  
 |  add_suffix(self, suffix)
 |      Suffix labels with string `suffix`.
 |      
 |      For Series, the row labels are suffixed.
 |      For DataFrame, the column labels are suffixed.
 |      
 |      Parameters
 |      ----------
 |      suffix : str
 |          The string to add after each label.
 |      
 |      Returns
 |      -------
 |      Series or DataFrame
 |          New Series or DataFrame with updated labels.
 |      
 |      See Also
 |      --------
 |      Series.add_prefix: Prefix row labels with string `prefix`.
 |      DataFrame.add_prefix: Prefix column labels with string `prefix`.
 |      
 |      Examples
 |      --------
 |      >>> s = pd.Series([1, 2, 3, 4])
 |      >>> s
 |      0    1
 |      1    2
 |      2    3
 |      3    4
 |      dtype: int64
 |      
 |      >>> s.add_suffix('_item')
 |      0_item    1
 |      1_item    2
 |      2_item    3
 |      3_item    4
 |      dtype: int64
 |      
 |      >>> df = pd.DataFrame({'A': [1, 2, 3, 4],  'B': [3, 4, 5, 6]})
 |      >>> df
 |         A  B
 |      0  1  3
 |      1  2  4
 |      2  3  5
 |      3  4  6
 |      
 |      >>> df.add_suffix('_col')
 |           A_col  B_col
 |      0       1       3
 |      1       2       4
 |      2       3       5
 |      3       4       6
 |  
 |  as_blocks(self, copy=True)
 |      Convert the frame to a dict of dtype -> Constructor Types that each has
 |      a homogeneous dtype.
 |      
 |      .. deprecated:: 0.21.0
 |      
 |      NOTE: the dtypes of the blocks WILL BE PRESERVED HERE (unlike in
 |            as_matrix)
 |      
 |      Parameters
 |      ----------
 |      copy : boolean, default True
 |      
 |      Returns
 |      -------
 |      values : a dict of dtype -> Constructor Types
 |  
 |  as_matrix(self, columns=None)
 |      Convert the frame to its Numpy-array representation.
 |      
 |      .. deprecated:: 0.23.0
 |          Use :meth:`DataFrame.values` instead.
 |      
 |      Parameters
 |      ----------
 |      columns: list, optional, default:None
 |          If None, return all columns, otherwise, returns specified columns.
 |      
 |      Returns
 |      -------
 |      values : ndarray
 |          If the caller is heterogeneous and contains booleans or objects,
 |          the result will be of dtype=object. See Notes.
 |      
 |      
 |      Notes
 |      -----
 |      Return is NOT a Numpy-matrix, rather, a Numpy-array.
 |      
 |      The dtype will be a lower-common-denominator dtype (implicit
 |      upcasting); that is to say if the dtypes (even of numeric types)
 |      are mixed, the one that accommodates all will be chosen. Use this
 |      with care if you are not dealing with the blocks.
 |      
 |      e.g. If the dtypes are float16 and float32, dtype will be upcast to
 |      float32.  If dtypes are int32 and uint8, dtype will be upcase to
 |      int32. By numpy.find_common_type convention, mixing int64 and uint64
 |      will result in a flot64 dtype.
 |      
 |      This method is provided for backwards compatibility. Generally,
 |      it is recommended to use '.values'.
 |      
 |      See Also
 |      --------
 |      pandas.DataFrame.values
 |  
 |  asfreq(self, freq, method=None, how=None, normalize=False, fill_value=None)
 |      Convert TimeSeries to specified frequency.
 |      
 |      Optionally provide filling method to pad/backfill missing values.
 |      
 |      Returns the original data conformed to a new index with the specified
 |      frequency. ``resample`` is more appropriate if an operation, such as
 |      summarization, is necessary to represent the data at the new frequency.
 |      
 |      Parameters
 |      ----------
 |      freq : DateOffset object, or string
 |      method : {'backfill'/'bfill', 'pad'/'ffill'}, default None
 |          Method to use for filling holes in reindexed Series (note this
 |          does not fill NaNs that already were present):
 |      
 |          * 'pad' / 'ffill': propagate last valid observation forward to next
 |            valid
 |          * 'backfill' / 'bfill': use NEXT valid observation to fill
 |      how : {'start', 'end'}, default end
 |          For PeriodIndex only, see PeriodIndex.asfreq
 |      normalize : bool, default False
 |          Whether to reset output index to midnight
 |      fill_value: scalar, optional
 |          Value to use for missing values, applied during upsampling (note
 |          this does not fill NaNs that already were present).
 |      
 |          .. versionadded:: 0.20.0
 |      
 |      Returns
 |      -------
 |      converted : type of caller
 |      
 |      Examples
 |      --------
 |      
 |      Start by creating a series with 4 one minute timestamps.
 |      
 |      >>> index = pd.date_range('1/1/2000', periods=4, freq='T')
 |      >>> series = pd.Series([0.0, None, 2.0, 3.0], index=index)
 |      >>> df = pd.DataFrame({'s':series})
 |      >>> df
 |                             s
 |      2000-01-01 00:00:00    0.0
 |      2000-01-01 00:01:00    NaN
 |      2000-01-01 00:02:00    2.0
 |      2000-01-01 00:03:00    3.0
 |      
 |      Upsample the series into 30 second bins.
 |      
 |      >>> df.asfreq(freq='30S')
 |                             s
 |      2000-01-01 00:00:00    0.0
 |      2000-01-01 00:00:30    NaN
 |      2000-01-01 00:01:00    NaN
 |      2000-01-01 00:01:30    NaN
 |      2000-01-01 00:02:00    2.0
 |      2000-01-01 00:02:30    NaN
 |      2000-01-01 00:03:00    3.0
 |      
 |      Upsample again, providing a ``fill value``.
 |      
 |      >>> df.asfreq(freq='30S', fill_value=9.0)
 |                             s
 |      2000-01-01 00:00:00    0.0
 |      2000-01-01 00:00:30    9.0
 |      2000-01-01 00:01:00    NaN
 |      2000-01-01 00:01:30    9.0
 |      2000-01-01 00:02:00    2.0
 |      2000-01-01 00:02:30    9.0
 |      2000-01-01 00:03:00    3.0
 |      
 |      Upsample again, providing a ``method``.
 |      
 |      >>> df.asfreq(freq='30S', method='bfill')
 |                             s
 |      2000-01-01 00:00:00    0.0
 |      2000-01-01 00:00:30    NaN
 |      2000-01-01 00:01:00    NaN
 |      2000-01-01 00:01:30    2.0
 |      2000-01-01 00:02:00    2.0
 |      2000-01-01 00:02:30    3.0
 |      2000-01-01 00:03:00    3.0
 |      
 |      See Also
 |      --------
 |      reindex
 |      
 |      Notes
 |      -----
 |      To learn more about the frequency strings, please see `this link
 |      <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__.
 |  
 |  asof(self, where, subset=None)
 |      The last row without any NaN is taken (or the last row without
 |      NaN considering only the subset of columns in the case of a DataFrame)
 |      
 |      .. versionadded:: 0.19.0 For DataFrame
 |      
 |      If there is no good value, NaN is returned for a Series
 |      a Series of NaN values for a DataFrame
 |      
 |      Parameters
 |      ----------
 |      where : date or array of dates
 |      subset : string or list of strings, default None
 |         if not None use these columns for NaN propagation
 |      
 |      Notes
 |      -----
 |      Dates are assumed to be sorted
 |      Raises if this is not the case
 |      
 |      Returns
 |      -------
 |      where is scalar
 |      
 |        - value or NaN if input is Series
 |        - Series if input is DataFrame
 |      
 |      where is Index: same shape object as input
 |      
 |      See Also
 |      --------
 |      merge_asof
 |  
 |  astype(self, dtype, copy=True, errors='raise', **kwargs)
 |      Cast a pandas object to a specified dtype ``dtype``.
 |      
 |      Parameters
 |      ----------
 |      dtype : data type, or dict of column name -> data type
 |          Use a numpy.dtype or Python type to cast entire pandas object to
 |          the same type. Alternatively, use {col: dtype, ...}, where col is a
 |          column label and dtype is a numpy.dtype or Python type to cast one
 |          or more of the DataFrame's columns to column-specific types.
 |      copy : bool, default True.
 |          Return a copy when ``copy=True`` (be very careful setting
 |          ``copy=False`` as changes to values then may propagate to other
 |          pandas objects).
 |      errors : {'raise', 'ignore'}, default 'raise'.
 |          Control raising of exceptions on invalid data for provided dtype.
 |      
 |          - ``raise`` : allow exceptions to be raised
 |          - ``ignore`` : suppress exceptions. On error return original object
 |      
 |          .. versionadded:: 0.20.0
 |      
 |      raise_on_error : raise on invalid input
 |          .. deprecated:: 0.20.0
 |             Use ``errors`` instead
 |      kwargs : keyword arguments to pass on to the constructor
 |      
 |      Returns
 |      -------
 |      casted : type of caller
 |      
 |      Examples
 |      --------
 |      >>> ser = pd.Series([1, 2], dtype='int32')
 |      >>> ser
 |      0    1
 |      1    2
 |      dtype: int32
 |      >>> ser.astype('int64')
 |      0    1
 |      1    2
 |      dtype: int64
 |      
 |      Convert to categorical type:
 |      
 |      >>> ser.astype('category')
 |      0    1
 |      1    2
 |      dtype: category
 |      Categories (2, int64): [1, 2]
 |      
 |      Convert to ordered categorical type with custom ordering:
 |      
 |      >>> ser.astype('category', ordered=True, categories=[2, 1])
 |      0    1
 |      1    2
 |      dtype: category
 |      Categories (2, int64): [2 < 1]
 |      
 |      Note that using ``copy=False`` and changing data on a new
 |      pandas object may propagate changes:
 |      
 |      >>> s1 = pd.Series([1,2])
 |      >>> s2 = s1.astype('int64', copy=False)
 |      >>> s2[0] = 10
 |      >>> s1  # note that s1[0] has changed too
 |      0    10
 |      1     2
 |      dtype: int64
 |      
 |      See also
 |      --------
 |      pandas.to_datetime : Convert argument to datetime.
 |      pandas.to_timedelta : Convert argument to timedelta.
 |      pandas.to_numeric : Convert argument to a numeric type.
 |      numpy.ndarray.astype : Cast a numpy array to a specified type.
 |  
 |  at_time(self, time, asof=False)
 |      Select values at particular time of day (e.g. 9:30AM).
 |      
 |      Raises
 |      ------
 |      TypeError
 |          If the index is not  a :class:`DatetimeIndex`
 |      
 |      Parameters
 |      ----------
 |      time : datetime.time or string
 |      
 |      Returns
 |      -------
 |      values_at_time : type of caller
 |      
 |      Examples
 |      --------
 |      >>> i = pd.date_range('2018-04-09', periods=4, freq='12H')
 |      >>> ts = pd.DataFrame({'A': [1,2,3,4]}, index=i)
 |      >>> ts
 |                           A
 |      2018-04-09 00:00:00  1
 |      2018-04-09 12:00:00  2
 |      2018-04-10 00:00:00  3
 |      2018-04-10 12:00:00  4
 |      
 |      >>> ts.at_time('12:00')
 |                           A
 |      2018-04-09 12:00:00  2
 |      2018-04-10 12:00:00  4
 |      
 |      See Also
 |      --------
 |      between_time : Select values between particular times of the day
 |      first : Select initial periods of time series based on a date offset
 |      last : Select final periods of time series based on a date offset
 |      DatetimeIndex.indexer_at_time : Get just the index locations for
 |          values at particular time of the day
 |  
 |  between_time(self, start_time, end_time, include_start=True, include_end=True)
 |      Select values between particular times of the day (e.g., 9:00-9:30 AM).
 |      
 |      By setting ``start_time`` to be later than ``end_time``,
 |      you can get the times that are *not* between the two times.
 |      
 |      Raises
 |      ------
 |      TypeError
 |          If the index is not  a :class:`DatetimeIndex`
 |      
 |      Parameters
 |      ----------
 |      start_time : datetime.time or string
 |      end_time : datetime.time or string
 |      include_start : boolean, default True
 |      include_end : boolean, default True
 |      
 |      Returns
 |      -------
 |      values_between_time : type of caller
 |      
 |      Examples
 |      --------
 |      >>> i = pd.date_range('2018-04-09', periods=4, freq='1D20min')
 |      >>> ts = pd.DataFrame({'A': [1,2,3,4]}, index=i)
 |      >>> ts
 |                           A
 |      2018-04-09 00:00:00  1
 |      2018-04-10 00:20:00  2
 |      2018-04-11 00:40:00  3
 |      2018-04-12 01:00:00  4
 |      
 |      >>> ts.between_time('0:15', '0:45')
 |                           A
 |      2018-04-10 00:20:00  2
 |      2018-04-11 00:40:00  3
 |      
 |      You get the times that are *not* between two times by setting
 |      ``start_time`` later than ``end_time``:
 |      
 |      >>> ts.between_time('0:45', '0:15')
 |                           A
 |      2018-04-09 00:00:00  1
 |      2018-04-12 01:00:00  4
 |      
 |      See Also
 |      --------
 |      at_time : Select values at a particular time of the day
 |      first : Select initial periods of time series based on a date offset
 |      last : Select final periods of time series based on a date offset
 |      DatetimeIndex.indexer_between_time : Get just the index locations for
 |          values between particular times of the day
 |  
 |  bfill(self, axis=None, inplace=False, limit=None, downcast=None)
 |      Synonym for :meth:`DataFrame.fillna(method='bfill') <DataFrame.fillna>`
 |  
 |  bool(self)
 |      Return the bool of a single element PandasObject.
 |      
 |      This must be a boolean scalar value, either True or False.  Raise a
 |      ValueError if the PandasObject does not have exactly 1 element, or that
 |      element is not boolean
 |  
 |  clip(self, lower=None, upper=None, axis=None, inplace=False, *args, **kwargs)
 |      Trim values at input threshold(s).
 |      
 |      Assigns values outside boundary to boundary values. Thresholds
 |      can be singular values or array like, and in the latter case
 |      the clipping is performed element-wise in the specified axis.
 |      
 |      Parameters
 |      ----------
 |      lower : float or array_like, default None
 |          Minimum threshold value. All values below this
 |          threshold will be set to it.
 |      upper : float or array_like, default None
 |          Maximum threshold value. All values above this
 |          threshold will be set to it.
 |      axis : int or string axis name, optional
 |          Align object with lower and upper along the given axis.
 |      inplace : boolean, default False
 |          Whether to perform the operation in place on the data.
 |      
 |          .. versionadded:: 0.21.0
 |      *args, **kwargs
 |          Additional keywords have no effect but might be accepted
 |          for compatibility with numpy.
 |      
 |      See Also
 |      --------
 |      clip_lower : Clip values below specified threshold(s).
 |      clip_upper : Clip values above specified threshold(s).
 |      
 |      Returns
 |      -------
 |      Series or DataFrame
 |          Same type as calling object with the values outside the
 |          clip boundaries replaced
 |      
 |      Examples
 |      --------
 |      >>> data = {'col_0': [9, -3, 0, -1, 5], 'col_1': [-2, -7, 6, 8, -5]}
 |      >>> df = pd.DataFrame(data)
 |      >>> df
 |         col_0  col_1
 |      0      9     -2
 |      1     -3     -7
 |      2      0      6
 |      3     -1      8
 |      4      5     -5
 |      
 |      Clips per column using lower and upper thresholds:
 |      
 |      >>> df.clip(-4, 6)
 |         col_0  col_1
 |      0      6     -2
 |      1     -3     -4
 |      2      0      6
 |      3     -1      6
 |      4      5     -4
 |      
 |      Clips using specific lower and upper thresholds per column element:
 |      
 |      >>> t = pd.Series([2, -4, -1, 6, 3])
 |      >>> t
 |      0    2
 |      1   -4
 |      2   -1
 |      3    6
 |      4    3
 |      dtype: int64
 |      
 |      >>> df.clip(t, t + 4, axis=0)
 |         col_0  col_1
 |      0      6      2
 |      1     -3     -4
 |      2      0      3
 |      3      6      8
 |      4      5      3
 |  
 |  clip_lower(self, threshold, axis=None, inplace=False)
 |      Return copy of the input with values below a threshold truncated.
 |      
 |      Parameters
 |      ----------
 |      threshold : numeric or array-like
 |          Minimum value allowed. All values below threshold will be set to
 |          this value.
 |      
 |          * float : every value is compared to `threshold`.
 |          * array-like : The shape of `threshold` should match the object
 |            it's compared to. When `self` is a Series, `threshold` should be
 |            the length. When `self` is a DataFrame, `threshold` should 2-D
 |            and the same shape as `self` for ``axis=None``, or 1-D and the
 |            same length as the axis being compared.
 |      
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          Align `self` with `threshold` along the given axis.
 |      
 |      inplace : boolean, default False
 |          Whether to perform the operation in place on the data.
 |      
 |          .. versionadded:: 0.21.0
 |      
 |      See Also
 |      --------
 |      Series.clip : Return copy of input with values below and above
 |          thresholds truncated.
 |      Series.clip_upper : Return copy of input with values above
 |          threshold truncated.
 |      
 |      Returns
 |      -------
 |      clipped : same type as input
 |      
 |      Examples
 |      --------
 |      Series single threshold clipping:
 |      
 |      >>> s = pd.Series([5, 6, 7, 8, 9])
 |      >>> s.clip_lower(8)
 |      0    8
 |      1    8
 |      2    8
 |      3    8
 |      4    9
 |      dtype: int64
 |      
 |      Series clipping element-wise using an array of thresholds. `threshold`
 |      should be the same length as the Series.
 |      
 |      >>> elemwise_thresholds = [4, 8, 7, 2, 5]
 |      >>> s.clip_lower(elemwise_thresholds)
 |      0    5
 |      1    8
 |      2    7
 |      3    8
 |      4    9
 |      dtype: int64
 |      
 |      DataFrames can be compared to a scalar.
 |      
 |      >>> df = pd.DataFrame({"A": [1, 3, 5], "B": [2, 4, 6]})
 |      >>> df
 |         A  B
 |      0  1  2
 |      1  3  4
 |      2  5  6
 |      
 |      >>> df.clip_lower(3)
 |         A  B
 |      0  3  3
 |      1  3  4
 |      2  5  6
 |      
 |      Or to an array of values. By default, `threshold` should be the same
 |      shape as the DataFrame.
 |      
 |      >>> df.clip_lower(np.array([[3, 4], [2, 2], [6, 2]]))
 |         A  B
 |      0  3  4
 |      1  3  4
 |      2  6  6
 |      
 |      Control how `threshold` is broadcast with `axis`. In this case
 |      `threshold` should be the same length as the axis specified by
 |      `axis`.
 |      
 |      >>> df.clip_lower(np.array([3, 3, 5]), axis='index')
 |         A  B
 |      0  3  3
 |      1  3  4
 |      2  5  6
 |      
 |      >>> df.clip_lower(np.array([4, 5]), axis='columns')
 |         A  B
 |      0  4  5
 |      1  4  5
 |      2  5  6
 |  
 |  clip_upper(self, threshold, axis=None, inplace=False)
 |      Return copy of input with values above given value(s) truncated.
 |      
 |      Parameters
 |      ----------
 |      threshold : float or array_like
 |      axis : int or string axis name, optional
 |          Align object with threshold along the given axis.
 |      inplace : boolean, default False
 |          Whether to perform the operation in place on the data
 |      
 |          .. versionadded:: 0.21.0
 |      
 |      See Also
 |      --------
 |      clip
 |      
 |      Returns
 |      -------
 |      clipped : same type as input
 |  
 |  consolidate(self, inplace=False)
 |      Compute NDFrame with "consolidated" internals (data of each dtype
 |      grouped together in a single ndarray).
 |      
 |      .. deprecated:: 0.20.0
 |          Consolidate will be an internal implementation only.
 |  
 |  convert_objects(self, convert_dates=True, convert_numeric=False, convert_timedeltas=True, copy=True)
 |      Attempt to infer better dtype for object columns.
 |      
 |      .. deprecated:: 0.21.0
 |      
 |      Parameters
 |      ----------
 |      convert_dates : boolean, default True
 |          If True, convert to date where possible. If 'coerce', force
 |          conversion, with unconvertible values becoming NaT.
 |      convert_numeric : boolean, default False
 |          If True, attempt to coerce to numbers (including strings), with
 |          unconvertible values becoming NaN.
 |      convert_timedeltas : boolean, default True
 |          If True, convert to timedelta where possible. If 'coerce', force
 |          conversion, with unconvertible values becoming NaT.
 |      copy : boolean, default True
 |          If True, return a copy even if no copy is necessary (e.g. no
 |          conversion was done). Note: This is meant for internal use, and
 |          should not be confused with inplace.
 |      
 |      See Also
 |      --------
 |      pandas.to_datetime : Convert argument to datetime.
 |      pandas.to_timedelta : Convert argument to timedelta.
 |      pandas.to_numeric : Return a fixed frequency timedelta index,
 |          with day as the default.
 |      
 |      Returns
 |      -------
 |      converted : same as input object
 |  
 |  copy(self, deep=True)
 |      Make a copy of this object's indices and data.
 |      
 |      When ``deep=True`` (default), a new object will be created with a
 |      copy of the calling object's data and indices. Modifications to
 |      the data or indices of the copy will not be reflected in the
 |      original object (see notes below).
 |      
 |      When ``deep=False``, a new object will be created without copying
 |      the calling object's data or index (only references to the data
 |      and index are copied). Any changes to the data of the original
 |      will be reflected in the shallow copy (and vice versa).
 |      
 |      Parameters
 |      ----------
 |      deep : bool, default True
 |          Make a deep copy, including a copy of the data and the indices.
 |          With ``deep=False`` neither the indices nor the data are copied.
 |      
 |      Returns
 |      -------
 |      copy : Series, DataFrame or Panel
 |          Object type matches caller.
 |      
 |      Notes
 |      -----
 |      When ``deep=True``, data is copied but actual Python objects
 |      will not be copied recursively, only the reference to the object.
 |      This is in contrast to `copy.deepcopy` in the Standard Library,
 |      which recursively copies object data (see examples below).
 |      
 |      While ``Index`` objects are copied when ``deep=True``, the underlying
 |      numpy array is not copied for performance reasons. Since ``Index`` is
 |      immutable, the underlying data can be safely shared and a copy
 |      is not needed.
 |      
 |      Examples
 |      --------
 |      >>> s = pd.Series([1, 2], index=["a", "b"])
 |      >>> s
 |      a    1
 |      b    2
 |      dtype: int64
 |      
 |      >>> s_copy = s.copy()
 |      >>> s_copy
 |      a    1
 |      b    2
 |      dtype: int64
 |      
 |      **Shallow copy versus default (deep) copy:**
 |      
 |      >>> s = pd.Series([1, 2], index=["a", "b"])
 |      >>> deep = s.copy()
 |      >>> shallow = s.copy(deep=False)
 |      
 |      Shallow copy shares data and index with original.
 |      
 |      >>> s is shallow
 |      False
 |      >>> s.values is shallow.values and s.index is shallow.index
 |      True
 |      
 |      Deep copy has own copy of data and index.
 |      
 |      >>> s is deep
 |      False
 |      >>> s.values is deep.values or s.index is deep.index
 |      False
 |      
 |      Updates to the data shared by shallow copy and original is reflected
 |      in both; deep copy remains unchanged.
 |      
 |      >>> s[0] = 3
 |      >>> shallow[1] = 4
 |      >>> s
 |      a    3
 |      b    4
 |      dtype: int64
 |      >>> shallow
 |      a    3
 |      b    4
 |      dtype: int64
 |      >>> deep
 |      a    1
 |      b    2
 |      dtype: int64
 |      
 |      Note that when copying an object containing Python objects, a deep copy
 |      will copy the data, but will not do so recursively. Updating a nested
 |      data object will be reflected in the deep copy.
 |      
 |      >>> s = pd.Series([[1, 2], [3, 4]])
 |      >>> deep = s.copy()
 |      >>> s[0][0] = 10
 |      >>> s
 |      0    [10, 2]
 |      1     [3, 4]
 |      dtype: object
 |      >>> deep
 |      0    [10, 2]
 |      1     [3, 4]
 |      dtype: object
 |  
 |  describe(self, percentiles=None, include=None, exclude=None)
 |      Generates descriptive statistics that summarize the central tendency,
 |      dispersion and shape of a dataset's distribution, excluding
 |      ``NaN`` values.
 |      
 |      Analyzes both numeric and object series, as well
 |      as ``DataFrame`` column sets of mixed data types. The output
 |      will vary depending on what is provided. Refer to the notes
 |      below for more detail.
 |      
 |      Parameters
 |      ----------
 |      percentiles : list-like of numbers, optional
 |          The percentiles to include in the output. All should
 |          fall between 0 and 1. The default is
 |          ``[.25, .5, .75]``, which returns the 25th, 50th, and
 |          75th percentiles.
 |      include : 'all', list-like of dtypes or None (default), optional
 |          A white list of data types to include in the result. Ignored
 |          for ``Series``. Here are the options:
 |      
 |          - 'all' : All columns of the input will be included in the output.
 |          - A list-like of dtypes : Limits the results to the
 |            provided data types.
 |            To limit the result to numeric types submit
 |            ``numpy.number``. To limit it instead to object columns submit
 |            the ``numpy.object`` data type. Strings
 |            can also be used in the style of
 |            ``select_dtypes`` (e.g. ``df.describe(include=['O'])``). To
 |            select pandas categorical columns, use ``'category'``
 |          - None (default) : The result will include all numeric columns.
 |      exclude : list-like of dtypes or None (default), optional,
 |          A black list of data types to omit from the result. Ignored
 |          for ``Series``. Here are the options:
 |      
 |          - A list-like of dtypes : Excludes the provided data types
 |            from the result. To exclude numeric types submit
 |            ``numpy.number``. To exclude object columns submit the data
 |            type ``numpy.object``. Strings can also be used in the style of
 |            ``select_dtypes`` (e.g. ``df.describe(include=['O'])``). To
 |            exclude pandas categorical columns, use ``'category'``
 |          - None (default) : The result will exclude nothing.
 |      
 |      Returns
 |      -------
 |      summary:  Series/DataFrame of summary statistics
 |      
 |      Notes
 |      -----
 |      For numeric data, the result's index will include ``count``,
 |      ``mean``, ``std``, ``min``, ``max`` as well as lower, ``50`` and
 |      upper percentiles. By default the lower percentile is ``25`` and the
 |      upper percentile is ``75``. The ``50`` percentile is the
 |      same as the median.
 |      
 |      For object data (e.g. strings or timestamps), the result's index
 |      will include ``count``, ``unique``, ``top``, and ``freq``. The ``top``
 |      is the most common value. The ``freq`` is the most common value's
 |      frequency. Timestamps also include the ``first`` and ``last`` items.
 |      
 |      If multiple object values have the highest count, then the
 |      ``count`` and ``top`` results will be arbitrarily chosen from
 |      among those with the highest count.
 |      
 |      For mixed data types provided via a ``DataFrame``, the default is to
 |      return only an analysis of numeric columns. If the dataframe consists
 |      only of object and categorical data without any numeric columns, the
 |      default is to return an analysis of both the object and categorical
 |      columns. If ``include='all'`` is provided as an option, the result
 |      will include a union of attributes of each type.
 |      
 |      The `include` and `exclude` parameters can be used to limit
 |      which columns in a ``DataFrame`` are analyzed for the output.
 |      The parameters are ignored when analyzing a ``Series``.
 |      
 |      Examples
 |      --------
 |      Describing a numeric ``Series``.
 |      
 |      >>> s = pd.Series([1, 2, 3])
 |      >>> s.describe()
 |      count    3.0
 |      mean     2.0
 |      std      1.0
 |      min      1.0
 |      25%      1.5
 |      50%      2.0
 |      75%      2.5
 |      max      3.0
 |      
 |      Describing a categorical ``Series``.
 |      
 |      >>> s = pd.Series(['a', 'a', 'b', 'c'])
 |      >>> s.describe()
 |      count     4
 |      unique    3
 |      top       a
 |      freq      2
 |      dtype: object
 |      
 |      Describing a timestamp ``Series``.
 |      
 |      >>> s = pd.Series([
 |      ...   np.datetime64("2000-01-01"),
 |      ...   np.datetime64("2010-01-01"),
 |      ...   np.datetime64("2010-01-01")
 |      ... ])
 |      >>> s.describe()
 |      count                       3
 |      unique                      2
 |      top       2010-01-01 00:00:00
 |      freq                        2
 |      first     2000-01-01 00:00:00
 |      last      2010-01-01 00:00:00
 |      dtype: object
 |      
 |      Describing a ``DataFrame``. By default only numeric fields
 |      are returned.
 |      
 |      >>> df = pd.DataFrame({ 'object': ['a', 'b', 'c'],
 |      ...                     'numeric': [1, 2, 3],
 |      ...                     'categorical': pd.Categorical(['d','e','f'])
 |      ...                   })
 |      >>> df.describe()
 |             numeric
 |      count      3.0
 |      mean       2.0
 |      std        1.0
 |      min        1.0
 |      25%        1.5
 |      50%        2.0
 |      75%        2.5
 |      max        3.0
 |      
 |      Describing all columns of a ``DataFrame`` regardless of data type.
 |      
 |      >>> df.describe(include='all')
 |              categorical  numeric object
 |      count            3      3.0      3
 |      unique           3      NaN      3
 |      top              f      NaN      c
 |      freq             1      NaN      1
 |      mean           NaN      2.0    NaN
 |      std            NaN      1.0    NaN
 |      min            NaN      1.0    NaN
 |      25%            NaN      1.5    NaN
 |      50%            NaN      2.0    NaN
 |      75%            NaN      2.5    NaN
 |      max            NaN      3.0    NaN
 |      
 |      Describing a column from a ``DataFrame`` by accessing it as
 |      an attribute.
 |      
 |      >>> df.numeric.describe()
 |      count    3.0
 |      mean     2.0
 |      std      1.0
 |      min      1.0
 |      25%      1.5
 |      50%      2.0
 |      75%      2.5
 |      max      3.0
 |      Name: numeric, dtype: float64
 |      
 |      Including only numeric columns in a ``DataFrame`` description.
 |      
 |      >>> df.describe(include=[np.number])
 |             numeric
 |      count      3.0
 |      mean       2.0
 |      std        1.0
 |      min        1.0
 |      25%        1.5
 |      50%        2.0
 |      75%        2.5
 |      max        3.0
 |      
 |      Including only string columns in a ``DataFrame`` description.
 |      
 |      >>> df.describe(include=[np.object])
 |             object
 |      count       3
 |      unique      3
 |      top         c
 |      freq        1
 |      
 |      Including only categorical columns from a ``DataFrame`` description.
 |      
 |      >>> df.describe(include=['category'])
 |             categorical
 |      count            3
 |      unique           3
 |      top              f
 |      freq             1
 |      
 |      Excluding numeric columns from a ``DataFrame`` description.
 |      
 |      >>> df.describe(exclude=[np.number])
 |             categorical object
 |      count            3      3
 |      unique           3      3
 |      top              f      c
 |      freq             1      1
 |      
 |      Excluding object columns from a ``DataFrame`` description.
 |      
 |      >>> df.describe(exclude=[np.object])
 |              categorical  numeric
 |      count            3      3.0
 |      unique           3      NaN
 |      top              f      NaN
 |      freq             1      NaN
 |      mean           NaN      2.0
 |      std            NaN      1.0
 |      min            NaN      1.0
 |      25%            NaN      1.5
 |      50%            NaN      2.0
 |      75%            NaN      2.5
 |      max            NaN      3.0
 |      
 |      See Also
 |      --------
 |      DataFrame.count
 |      DataFrame.max
 |      DataFrame.min
 |      DataFrame.mean
 |      DataFrame.std
 |      DataFrame.select_dtypes
 |  
 |  equals(self, other)
 |      Determines if two NDFrame objects contain the same elements. NaNs in
 |      the same location are considered equal.
 |  
 |  ffill(self, axis=None, inplace=False, limit=None, downcast=None)
 |      Synonym for :meth:`DataFrame.fillna(method='ffill') <DataFrame.fillna>`
 |  
 |  filter(self, items=None, like=None, regex=None, axis=None)
 |      Subset rows or columns of dataframe according to labels in
 |      the specified index.
 |      
 |      Note that this routine does not filter a dataframe on its
 |      contents. The filter is applied to the labels of the index.
 |      
 |      Parameters
 |      ----------
 |      items : list-like
 |          List of info axis to restrict to (must not all be present)
 |      like : string
 |          Keep info axis where "arg in col == True"
 |      regex : string (regular expression)
 |          Keep info axis with re.search(regex, col) == True
 |      axis : int or string axis name
 |          The axis to filter on.  By default this is the info axis,
 |          'index' for Series, 'columns' for DataFrame
 |      
 |      Returns
 |      -------
 |      same type as input object
 |      
 |      Examples
 |      --------
 |      >>> df
 |      one  two  three
 |      mouse     1    2      3
 |      rabbit    4    5      6
 |      
 |      >>> # select columns by name
 |      >>> df.filter(items=['one', 'three'])
 |      one  three
 |      mouse     1      3
 |      rabbit    4      6
 |      
 |      >>> # select columns by regular expression
 |      >>> df.filter(regex='e$', axis=1)
 |      one  three
 |      mouse     1      3
 |      rabbit    4      6
 |      
 |      >>> # select rows containing 'bbi'
 |      >>> df.filter(like='bbi', axis=0)
 |      one  two  three
 |      rabbit    4    5      6
 |      
 |      See Also
 |      --------
 |      pandas.DataFrame.loc
 |      
 |      Notes
 |      -----
 |      The ``items``, ``like``, and ``regex`` parameters are
 |      enforced to be mutually exclusive.
 |      
 |      ``axis`` defaults to the info axis that is used when indexing
 |      with ``[]``.
 |  
 |  first(self, offset)
 |      Convenience method for subsetting initial periods of time series data
 |      based on a date offset.
 |      
 |      Raises
 |      ------
 |      TypeError
 |          If the index is not  a :class:`DatetimeIndex`
 |      
 |      Parameters
 |      ----------
 |      offset : string, DateOffset, dateutil.relativedelta
 |      
 |      Examples
 |      --------
 |      >>> i = pd.date_range('2018-04-09', periods=4, freq='2D')
 |      >>> ts = pd.DataFrame({'A': [1,2,3,4]}, index=i)
 |      >>> ts
 |                  A
 |      2018-04-09  1
 |      2018-04-11  2
 |      2018-04-13  3
 |      2018-04-15  4
 |      
 |      Get the rows for the first 3 days:
 |      
 |      >>> ts.first('3D')
 |                  A
 |      2018-04-09  1
 |      2018-04-11  2
 |      
 |      Notice the data for 3 first calender days were returned, not the first
 |      3 days observed in the dataset, and therefore data for 2018-04-13 was
 |      not returned.
 |      
 |      Returns
 |      -------
 |      subset : type of caller
 |      
 |      See Also
 |      --------
 |      last : Select final periods of time series based on a date offset
 |      at_time : Select values at a particular time of the day
 |      between_time : Select values between particular times of the day
 |  
 |  first_valid_index(self)
 |      Return index for first non-NA/null value.
 |      
 |      Notes
 |      --------
 |      If all elements are non-NA/null, returns None.
 |      Also returns None for empty NDFrame.
 |      
 |      Returns
 |      --------
 |      scalar : type of index
 |  
 |  get(self, key, default=None)
 |      Get item from object for given key (DataFrame column, Panel slice,
 |      etc.). Returns default value if not found.
 |      
 |      Parameters
 |      ----------
 |      key : object
 |      
 |      Returns
 |      -------
 |      value : type of items contained in object
 |  
 |  get_dtype_counts(self)
 |      Return counts of unique dtypes in this object.
 |      
 |      Returns
 |      -------
 |      dtype : Series
 |          Series with the count of columns with each dtype.
 |      
 |      See Also
 |      --------
 |      dtypes : Return the dtypes in this object.
 |      
 |      Examples
 |      --------
 |      >>> a = [['a', 1, 1.0], ['b', 2, 2.0], ['c', 3, 3.0]]
 |      >>> df = pd.DataFrame(a, columns=['str', 'int', 'float'])
 |      >>> df
 |        str  int  float
 |      0   a    1    1.0
 |      1   b    2    2.0
 |      2   c    3    3.0
 |      
 |      >>> df.get_dtype_counts()
 |      float64    1
 |      int64      1
 |      object     1
 |      dtype: int64
 |  
 |  get_ftype_counts(self)
 |      Return counts of unique ftypes in this object.
 |      
 |      .. deprecated:: 0.23.0
 |      
 |      This is useful for SparseDataFrame or for DataFrames containing
 |      sparse arrays.
 |      
 |      Returns
 |      -------
 |      dtype : Series
 |          Series with the count of columns with each type and
 |          sparsity (dense/sparse)
 |      
 |      See Also
 |      --------
 |      ftypes : Return ftypes (indication of sparse/dense and dtype) in
 |          this object.
 |      
 |      Examples
 |      --------
 |      >>> a = [['a', 1, 1.0], ['b', 2, 2.0], ['c', 3, 3.0]]
 |      >>> df = pd.DataFrame(a, columns=['str', 'int', 'float'])
 |      >>> df
 |        str  int  float
 |      0   a    1    1.0
 |      1   b    2    2.0
 |      2   c    3    3.0
 |      
 |      >>> df.get_ftype_counts()
 |      float64:dense    1
 |      int64:dense      1
 |      object:dense     1
 |      dtype: int64
 |  
 |  get_values(self)
 |      Return an ndarray after converting sparse values to dense.
 |      
 |      This is the same as ``.values`` for non-sparse data. For sparse
 |      data contained in a `pandas.SparseArray`, the data are first
 |      converted to a dense representation.
 |      
 |      Returns
 |      -------
 |      numpy.ndarray
 |          Numpy representation of DataFrame
 |      
 |      See Also
 |      --------
 |      values : Numpy representation of DataFrame.
 |      pandas.SparseArray : Container for sparse data.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'a': [1, 2], 'b': [True, False],
 |      ...                    'c': [1.0, 2.0]})
 |      >>> df
 |         a      b    c
 |      0  1   True  1.0
 |      1  2  False  2.0
 |      
 |      >>> df.get_values()
 |      array([[1, True, 1.0], [2, False, 2.0]], dtype=object)
 |      
 |      >>> df = pd.DataFrame({"a": pd.SparseArray([1, None, None]),
 |      ...                    "c": [1.0, 2.0, 3.0]})
 |      >>> df
 |           a    c
 |      0  1.0  1.0
 |      1  NaN  2.0
 |      2  NaN  3.0
 |      
 |      >>> df.get_values()
 |      array([[ 1.,  1.],
 |             [nan,  2.],
 |             [nan,  3.]])
 |  
 |  groupby(self, by=None, axis=0, level=None, as_index=True, sort=True, group_keys=True, squeeze=False, observed=False, **kwargs)
 |      Group series using mapper (dict or key function, apply given function
 |      to group, return result as series) or by a series of columns.
 |      
 |      Parameters
 |      ----------
 |      by : mapping, function, label, or list of labels
 |          Used to determine the groups for the groupby.
 |          If ``by`` is a function, it's called on each value of the object's
 |          index. If a dict or Series is passed, the Series or dict VALUES
 |          will be used to determine the groups (the Series' values are first
 |          aligned; see ``.align()`` method). If an ndarray is passed, the
 |          values are used as-is determine the groups. A label or list of
 |          labels may be passed to group by the columns in ``self``. Notice
 |          that a tuple is interpreted a (single) key.
 |      axis : int, default 0
 |      level : int, level name, or sequence of such, default None
 |          If the axis is a MultiIndex (hierarchical), group by a particular
 |          level or levels
 |      as_index : boolean, default True
 |          For aggregated output, return object with group labels as the
 |          index. Only relevant for DataFrame input. as_index=False is
 |          effectively "SQL-style" grouped output
 |      sort : boolean, default True
 |          Sort group keys. Get better performance by turning this off.
 |          Note this does not influence the order of observations within each
 |          group.  groupby preserves the order of rows within each group.
 |      group_keys : boolean, default True
 |          When calling apply, add group keys to index to identify pieces
 |      squeeze : boolean, default False
 |          reduce the dimensionality of the return type if possible,
 |          otherwise return a consistent type
 |      observed : boolean, default False
 |          This only applies if any of the groupers are Categoricals
 |          If True: only show observed values for categorical groupers.
 |          If False: show all values for categorical groupers.
 |      
 |          .. versionadded:: 0.23.0
 |      
 |      Returns
 |      -------
 |      GroupBy object
 |      
 |      Examples
 |      --------
 |      DataFrame results
 |      
 |      >>> data.groupby(func, axis=0).mean()
 |      >>> data.groupby(['col1', 'col2'])['col3'].mean()
 |      
 |      DataFrame with hierarchical index
 |      
 |      >>> data.groupby(['col1', 'col2']).mean()
 |      
 |      Notes
 |      -----
 |      See the `user guide
 |      <http://pandas.pydata.org/pandas-docs/stable/groupby.html>`_ for more.
 |      
 |      See also
 |      --------
 |      resample : Convenience method for frequency conversion and resampling
 |          of time series.
 |  
 |  head(self, n=5)
 |      Return the first `n` rows.
 |      
 |      This function returns the first `n` rows for the object based
 |      on position. It is useful for quickly testing if your object
 |      has the right type of data in it.
 |      
 |      Parameters
 |      ----------
 |      n : int, default 5
 |          Number of rows to select.
 |      
 |      Returns
 |      -------
 |      obj_head : type of caller
 |          The first `n` rows of the caller object.
 |      
 |      See Also
 |      --------
 |      pandas.DataFrame.tail: Returns the last `n` rows.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'animal':['alligator', 'bee', 'falcon', 'lion',
 |      ...                    'monkey', 'parrot', 'shark', 'whale', 'zebra']})
 |      >>> df
 |            animal
 |      0  alligator
 |      1        bee
 |      2     falcon
 |      3       lion
 |      4     monkey
 |      5     parrot
 |      6      shark
 |      7      whale
 |      8      zebra
 |      
 |      Viewing the first 5 lines
 |      
 |      >>> df.head()
 |            animal
 |      0  alligator
 |      1        bee
 |      2     falcon
 |      3       lion
 |      4     monkey
 |      
 |      Viewing the first `n` lines (three in this case)
 |      
 |      >>> df.head(3)
 |            animal
 |      0  alligator
 |      1        bee
 |      2     falcon
 |  
 |  infer_objects(self)
 |      Attempt to infer better dtypes for object columns.
 |      
 |      Attempts soft conversion of object-dtyped
 |      columns, leaving non-object and unconvertible
 |      columns unchanged. The inference rules are the
 |      same as during normal Series/DataFrame construction.
 |      
 |      .. versionadded:: 0.21.0
 |      
 |      See Also
 |      --------
 |      pandas.to_datetime : Convert argument to datetime.
 |      pandas.to_timedelta : Convert argument to timedelta.
 |      pandas.to_numeric : Convert argument to numeric typeR
 |      
 |      Returns
 |      -------
 |      converted : same type as input object
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({"A": ["a", 1, 2, 3]})
 |      >>> df = df.iloc[1:]
 |      >>> df
 |         A
 |      1  1
 |      2  2
 |      3  3
 |      
 |      >>> df.dtypes
 |      A    object
 |      dtype: object
 |      
 |      >>> df.infer_objects().dtypes
 |      A    int64
 |      dtype: object
 |  
 |  interpolate(self, method='linear', axis=0, limit=None, inplace=False, limit_direction='forward', limit_area=None, downcast=None, **kwargs)
 |      Interpolate values according to different methods.
 |      
 |      Please note that only ``method='linear'`` is supported for
 |      DataFrames/Series with a MultiIndex.
 |      
 |      Parameters
 |      ----------
 |      method : {'linear', 'time', 'index', 'values', 'nearest', 'zero',
 |                'slinear', 'quadratic', 'cubic', 'barycentric', 'krogh',
 |                'polynomial', 'spline', 'piecewise_polynomial',
 |                'from_derivatives', 'pchip', 'akima'}
 |      
 |          * 'linear': ignore the index and treat the values as equally
 |            spaced. This is the only method supported on MultiIndexes.
 |            default
 |          * 'time': interpolation works on daily and higher resolution
 |            data to interpolate given length of interval
 |          * 'index', 'values': use the actual numerical values of the index
 |          * 'nearest', 'zero', 'slinear', 'quadratic', 'cubic',
 |            'barycentric', 'polynomial' is passed to
 |            ``scipy.interpolate.interp1d``. Both 'polynomial' and 'spline'
 |            require that you also specify an `order` (int),
 |            e.g. df.interpolate(method='polynomial', order=4).
 |            These use the actual numerical values of the index.
 |          * 'krogh', 'piecewise_polynomial', 'spline', 'pchip' and 'akima'
 |            are all wrappers around the scipy interpolation methods of
 |            similar names. These use the actual numerical values of the
 |            index. For more information on their behavior, see the
 |            `scipy documentation
 |            <http://docs.scipy.org/doc/scipy/reference/interpolate.html#univariate-interpolation>`__
 |            and `tutorial documentation
 |            <http://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html>`__
 |          * 'from_derivatives' refers to BPoly.from_derivatives which
 |            replaces 'piecewise_polynomial' interpolation method in
 |            scipy 0.18
 |      
 |          .. versionadded:: 0.18.1
 |      
 |             Added support for the 'akima' method
 |             Added interpolate method 'from_derivatives' which replaces
 |             'piecewise_polynomial' in scipy 0.18; backwards-compatible with
 |             scipy < 0.18
 |      
 |      axis : {0, 1}, default 0
 |          * 0: fill column-by-column
 |          * 1: fill row-by-row
 |      limit : int, default None.
 |          Maximum number of consecutive NaNs to fill. Must be greater than 0.
 |      limit_direction : {'forward', 'backward', 'both'}, default 'forward'
 |      limit_area : {'inside', 'outside'}, default None
 |          * None: (default) no fill restriction
 |          * 'inside' Only fill NaNs surrounded by valid values (interpolate).
 |          * 'outside' Only fill NaNs outside valid values (extrapolate).
 |      
 |          If limit is specified, consecutive NaNs will be filled in this
 |          direction.
 |      
 |          .. versionadded:: 0.21.0
 |      inplace : bool, default False
 |          Update the NDFrame in place if possible.
 |      downcast : optional, 'infer' or None, defaults to None
 |          Downcast dtypes if possible.
 |      kwargs : keyword arguments to pass on to the interpolating function.
 |      
 |      Returns
 |      -------
 |      Series or DataFrame of same shape interpolated at the NaNs
 |      
 |      See Also
 |      --------
 |      reindex, replace, fillna
 |      
 |      Examples
 |      --------
 |      
 |      Filling in NaNs
 |      
 |      >>> s = pd.Series([0, 1, np.nan, 3])
 |      >>> s.interpolate()
 |      0    0
 |      1    1
 |      2    2
 |      3    3
 |      dtype: float64
 |  
 |  keys(self)
 |      Get the 'info axis' (see Indexing for more)
 |      
 |      This is index for Series, columns for DataFrame and major_axis for
 |      Panel.
 |  
 |  last(self, offset)
 |      Convenience method for subsetting final periods of time series data
 |      based on a date offset.
 |      
 |      Raises
 |      ------
 |      TypeError
 |          If the index is not  a :class:`DatetimeIndex`
 |      
 |      Parameters
 |      ----------
 |      offset : string, DateOffset, dateutil.relativedelta
 |      
 |      Examples
 |      --------
 |      >>> i = pd.date_range('2018-04-09', periods=4, freq='2D')
 |      >>> ts = pd.DataFrame({'A': [1,2,3,4]}, index=i)
 |      >>> ts
 |                  A
 |      2018-04-09  1
 |      2018-04-11  2
 |      2018-04-13  3
 |      2018-04-15  4
 |      
 |      Get the rows for the last 3 days:
 |      
 |      >>> ts.last('3D')
 |                  A
 |      2018-04-13  3
 |      2018-04-15  4
 |      
 |      Notice the data for 3 last calender days were returned, not the last
 |      3 observed days in the dataset, and therefore data for 2018-04-11 was
 |      not returned.
 |      
 |      Returns
 |      -------
 |      subset : type of caller
 |      
 |      See Also
 |      --------
 |      first : Select initial periods of time series based on a date offset
 |      at_time : Select values at a particular time of the day
 |      between_time : Select values between particular times of the day
 |  
 |  last_valid_index(self)
 |      Return index for last non-NA/null value.
 |      
 |      Notes
 |      --------
 |      If all elements are non-NA/null, returns None.
 |      Also returns None for empty NDFrame.
 |      
 |      Returns
 |      --------
 |      scalar : type of index
 |  
 |  mask(self, cond, other=nan, inplace=False, axis=None, level=None, errors='raise', try_cast=False, raise_on_error=None)
 |      Return an object of same shape as self and whose corresponding
 |      entries are from self where `cond` is False and otherwise are from
 |      `other`.
 |      
 |      Parameters
 |      ----------
 |      cond : boolean NDFrame, array-like, or callable
 |          Where `cond` is False, keep the original value. Where
 |          True, replace with corresponding value from `other`.
 |          If `cond` is callable, it is computed on the NDFrame and
 |          should return boolean NDFrame or array. The callable must
 |          not change input NDFrame (though pandas doesn't check it).
 |      
 |          .. versionadded:: 0.18.1
 |              A callable can be used as cond.
 |      
 |      other : scalar, NDFrame, or callable
 |          Entries where `cond` is True are replaced with
 |          corresponding value from `other`.
 |          If other is callable, it is computed on the NDFrame and
 |          should return scalar or NDFrame. The callable must not
 |          change input NDFrame (though pandas doesn't check it).
 |      
 |          .. versionadded:: 0.18.1
 |              A callable can be used as other.
 |      
 |      inplace : boolean, default False
 |          Whether to perform the operation in place on the data
 |      axis : alignment axis if needed, default None
 |      level : alignment level if needed, default None
 |      errors : str, {'raise', 'ignore'}, default 'raise'
 |          - ``raise`` : allow exceptions to be raised
 |          - ``ignore`` : suppress exceptions. On error return original object
 |      
 |          Note that currently this parameter won't affect
 |          the results and will always coerce to a suitable dtype.
 |      
 |      try_cast : boolean, default False
 |          try to cast the result back to the input type (if possible),
 |      raise_on_error : boolean, default True
 |          Whether to raise on invalid data types (e.g. trying to where on
 |          strings)
 |      
 |          .. deprecated:: 0.21.0
 |      
 |      Returns
 |      -------
 |      wh : same type as caller
 |      
 |      Notes
 |      -----
 |      The mask method is an application of the if-then idiom. For each
 |      element in the calling DataFrame, if ``cond`` is ``False`` the
 |      element is used; otherwise the corresponding element from the DataFrame
 |      ``other`` is used.
 |      
 |      The signature for :func:`DataFrame.where` differs from
 |      :func:`numpy.where`. Roughly ``df1.where(m, df2)`` is equivalent to
 |      ``np.where(m, df1, df2)``.
 |      
 |      For further details and examples see the ``mask`` documentation in
 |      :ref:`indexing <indexing.where_mask>`.
 |      
 |      Examples
 |      --------
 |      >>> s = pd.Series(range(5))
 |      >>> s.where(s > 0)
 |      0    NaN
 |      1    1.0
 |      2    2.0
 |      3    3.0
 |      4    4.0
 |      
 |      >>> s.mask(s > 0)
 |      0    0.0
 |      1    NaN
 |      2    NaN
 |      3    NaN
 |      4    NaN
 |      
 |      >>> s.where(s > 1, 10)
 |      0    10.0
 |      1    10.0
 |      2    2.0
 |      3    3.0
 |      4    4.0
 |      
 |      >>> df = pd.DataFrame(np.arange(10).reshape(-1, 2), columns=['A', 'B'])
 |      >>> m = df % 3 == 0
 |      >>> df.where(m, -df)
 |         A  B
 |      0  0 -1
 |      1 -2  3
 |      2 -4 -5
 |      3  6 -7
 |      4 -8  9
 |      >>> df.where(m, -df) == np.where(m, df, -df)
 |            A     B
 |      0  True  True
 |      1  True  True
 |      2  True  True
 |      3  True  True
 |      4  True  True
 |      >>> df.where(m, -df) == df.mask(~m, -df)
 |            A     B
 |      0  True  True
 |      1  True  True
 |      2  True  True
 |      3  True  True
 |      4  True  True
 |      
 |      See Also
 |      --------
 |      :func:`DataFrame.where`
 |  
 |  pct_change(self, periods=1, fill_method='pad', limit=None, freq=None, **kwargs)
 |      Percentage change between the current and a prior element.
 |      
 |      Computes the percentage change from the immediately previous row by
 |      default. This is useful in comparing the percentage of change in a time
 |      series of elements.
 |      
 |      Parameters
 |      ----------
 |      periods : int, default 1
 |          Periods to shift for forming percent change.
 |      fill_method : str, default 'pad'
 |          How to handle NAs before computing percent changes.
 |      limit : int, default None
 |          The number of consecutive NAs to fill before stopping.
 |      freq : DateOffset, timedelta, or offset alias string, optional
 |          Increment to use from time series API (e.g. 'M' or BDay()).
 |      **kwargs
 |          Additional keyword arguments are passed into
 |          `DataFrame.shift` or `Series.shift`.
 |      
 |      Returns
 |      -------
 |      chg : Series or DataFrame
 |          The same type as the calling object.
 |      
 |      See Also
 |      --------
 |      Series.diff : Compute the difference of two elements in a Series.
 |      DataFrame.diff : Compute the difference of two elements in a DataFrame.
 |      Series.shift : Shift the index by some number of periods.
 |      DataFrame.shift : Shift the index by some number of periods.
 |      
 |      Examples
 |      --------
 |      **Series**
 |      
 |      >>> s = pd.Series([90, 91, 85])
 |      >>> s
 |      0    90
 |      1    91
 |      2    85
 |      dtype: int64
 |      
 |      >>> s.pct_change()
 |      0         NaN
 |      1    0.011111
 |      2   -0.065934
 |      dtype: float64
 |      
 |      >>> s.pct_change(periods=2)
 |      0         NaN
 |      1         NaN
 |      2   -0.055556
 |      dtype: float64
 |      
 |      See the percentage change in a Series where filling NAs with last
 |      valid observation forward to next valid.
 |      
 |      >>> s = pd.Series([90, 91, None, 85])
 |      >>> s
 |      0    90.0
 |      1    91.0
 |      2     NaN
 |      3    85.0
 |      dtype: float64
 |      
 |      >>> s.pct_change(fill_method='ffill')
 |      0         NaN
 |      1    0.011111
 |      2    0.000000
 |      3   -0.065934
 |      dtype: float64
 |      
 |      **DataFrame**
 |      
 |      Percentage change in French franc, Deutsche Mark, and Italian lira from
 |      1980-01-01 to 1980-03-01.
 |      
 |      >>> df = pd.DataFrame({
 |      ...     'FR': [4.0405, 4.0963, 4.3149],
 |      ...     'GR': [1.7246, 1.7482, 1.8519],
 |      ...     'IT': [804.74, 810.01, 860.13]},
 |      ...     index=['1980-01-01', '1980-02-01', '1980-03-01'])
 |      >>> df
 |                      FR      GR      IT
 |      1980-01-01  4.0405  1.7246  804.74
 |      1980-02-01  4.0963  1.7482  810.01
 |      1980-03-01  4.3149  1.8519  860.13
 |      
 |      >>> df.pct_change()
 |                        FR        GR        IT
 |      1980-01-01       NaN       NaN       NaN
 |      1980-02-01  0.013810  0.013684  0.006549
 |      1980-03-01  0.053365  0.059318  0.061876
 |      
 |      Percentage of change in GOOG and APPL stock volume. Shows computing
 |      the percentage change between columns.
 |      
 |      >>> df = pd.DataFrame({
 |      ...     '2016': [1769950, 30586265],
 |      ...     '2015': [1500923, 40912316],
 |      ...     '2014': [1371819, 41403351]},
 |      ...     index=['GOOG', 'APPL'])
 |      >>> df
 |                2016      2015      2014
 |      GOOG   1769950   1500923   1371819
 |      APPL  30586265  40912316  41403351
 |      
 |      >>> df.pct_change(axis='columns')
 |            2016      2015      2014
 |      GOOG   NaN -0.151997 -0.086016
 |      APPL   NaN  0.337604  0.012002
 |  
 |  pipe(self, func, *args, **kwargs)
 |      Apply func(self, \*args, \*\*kwargs)
 |      
 |      Parameters
 |      ----------
 |      func : function
 |          function to apply to the NDFrame.
 |          ``args``, and ``kwargs`` are passed into ``func``.
 |          Alternatively a ``(callable, data_keyword)`` tuple where
 |          ``data_keyword`` is a string indicating the keyword of
 |          ``callable`` that expects the NDFrame.
 |      args : iterable, optional
 |          positional arguments passed into ``func``.
 |      kwargs : mapping, optional
 |          a dictionary of keyword arguments passed into ``func``.
 |      
 |      Returns
 |      -------
 |      object : the return type of ``func``.
 |      
 |      Notes
 |      -----
 |      
 |      Use ``.pipe`` when chaining together functions that expect
 |      Series, DataFrames or GroupBy objects. Instead of writing
 |      
 |      >>> f(g(h(df), arg1=a), arg2=b, arg3=c)
 |      
 |      You can write
 |      
 |      >>> (df.pipe(h)
 |      ...    .pipe(g, arg1=a)
 |      ...    .pipe(f, arg2=b, arg3=c)
 |      ... )
 |      
 |      If you have a function that takes the data as (say) the second
 |      argument, pass a tuple indicating which keyword expects the
 |      data. For example, suppose ``f`` takes its data as ``arg2``:
 |      
 |      >>> (df.pipe(h)
 |      ...    .pipe(g, arg1=a)
 |      ...    .pipe((f, 'arg2'), arg1=a, arg3=c)
 |      ...  )
 |      
 |      See Also
 |      --------
 |      pandas.DataFrame.apply
 |      pandas.DataFrame.applymap
 |      pandas.Series.map
 |  
 |  pop(self, item)
 |      Return item and drop from frame. Raise KeyError if not found.
 |      
 |      Parameters
 |      ----------
 |      item : str
 |          Column label to be popped
 |      
 |      Returns
 |      -------
 |      popped : Series
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame([('falcon', 'bird',    389.0),
 |      ...                    ('parrot', 'bird',     24.0),
 |      ...                    ('lion',   'mammal',   80.5),
 |      ...                    ('monkey', 'mammal', np.nan)],
 |      ...                   columns=('name', 'class', 'max_speed'))
 |      >>> df
 |           name   class  max_speed
 |      0  falcon    bird      389.0
 |      1  parrot    bird       24.0
 |      2    lion  mammal       80.5
 |      3  monkey  mammal        NaN
 |      
 |      >>> df.pop('class')
 |      0      bird
 |      1      bird
 |      2    mammal
 |      3    mammal
 |      Name: class, dtype: object
 |      
 |      >>> df
 |           name  max_speed
 |      0  falcon      389.0
 |      1  parrot       24.0
 |      2    lion       80.5
 |      3  monkey        NaN
 |  
 |  rank(self, axis=0, method='average', numeric_only=None, na_option='keep', ascending=True, pct=False)
 |      Compute numerical data ranks (1 through n) along axis. Equal values are
 |      assigned a rank that is the average of the ranks of those values
 |      
 |      Parameters
 |      ----------
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          index to direct ranking
 |      method : {'average', 'min', 'max', 'first', 'dense'}
 |          * average: average rank of group
 |          * min: lowest rank in group
 |          * max: highest rank in group
 |          * first: ranks assigned in order they appear in the array
 |          * dense: like 'min', but rank always increases by 1 between groups
 |      numeric_only : boolean, default None
 |          Include only float, int, boolean data. Valid only for DataFrame or
 |          Panel objects
 |      na_option : {'keep', 'top', 'bottom'}
 |          * keep: leave NA values where they are
 |          * top: smallest rank if ascending
 |          * bottom: smallest rank if descending
 |      ascending : boolean, default True
 |          False for ranks by high (1) to low (N)
 |      pct : boolean, default False
 |          Computes percentage rank of data
 |      
 |      Returns
 |      -------
 |      ranks : same type as caller
 |  
 |  reindex_like(self, other, method=None, copy=True, limit=None, tolerance=None)
 |      Return an object with matching indices to myself.
 |      
 |      Parameters
 |      ----------
 |      other : Object
 |      method : string or None
 |      copy : boolean, default True
 |      limit : int, default None
 |          Maximum number of consecutive labels to fill for inexact matches.
 |      tolerance : optional
 |          Maximum distance between labels of the other object and this
 |          object for inexact matches. Can be list-like.
 |      
 |          .. versionadded:: 0.21.0 (list-like tolerance)
 |      
 |      Notes
 |      -----
 |      Like calling s.reindex(index=other.index, columns=other.columns,
 |                             method=...)
 |      
 |      Returns
 |      -------
 |      reindexed : same as input
 |  
 |  rename_axis(self, mapper, axis=0, copy=True, inplace=False)
 |      Alter the name of the index or columns.
 |      
 |      Parameters
 |      ----------
 |      mapper : scalar, list-like, optional
 |          Value to set as the axis name attribute.
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          The index or the name of the axis.
 |      copy : boolean, default True
 |          Also copy underlying data.
 |      inplace : boolean, default False
 |          Modifies the object directly, instead of creating a new Series
 |          or DataFrame.
 |      
 |      Returns
 |      -------
 |      renamed : Series, DataFrame, or None
 |          The same type as the caller or None if `inplace` is True.
 |      
 |      Notes
 |      -----
 |      Prior to version 0.21.0, ``rename_axis`` could also be used to change
 |      the axis *labels* by passing a mapping or scalar. This behavior is
 |      deprecated and will be removed in a future version. Use ``rename``
 |      instead.
 |      
 |      See Also
 |      --------
 |      pandas.Series.rename : Alter Series index labels or name
 |      pandas.DataFrame.rename : Alter DataFrame index labels or name
 |      pandas.Index.rename : Set new names on index
 |      
 |      Examples
 |      --------
 |      **Series**
 |      
 |      >>> s = pd.Series([1, 2, 3])
 |      >>> s.rename_axis("foo")
 |      foo
 |      0    1
 |      1    2
 |      2    3
 |      dtype: int64
 |      
 |      **DataFrame**
 |      
 |      >>> df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
 |      >>> df.rename_axis("foo")
 |           A  B
 |      foo
 |      0    1  4
 |      1    2  5
 |      2    3  6
 |      
 |      >>> df.rename_axis("bar", axis="columns")
 |      bar  A  B
 |      0    1  4
 |      1    2  5
 |      2    3  6
 |  
 |  resample(self, rule, how=None, axis=0, fill_method=None, closed=None, label=None, convention='start', kind=None, loffset=None, limit=None, base=0, on=None, level=None)
 |      Convenience method for frequency conversion and resampling of time
 |      series.  Object must have a datetime-like index (DatetimeIndex,
 |      PeriodIndex, or TimedeltaIndex), or pass datetime-like values
 |      to the on or level keyword.
 |      
 |      Parameters
 |      ----------
 |      rule : string
 |          the offset string or object representing target conversion
 |      axis : int, optional, default 0
 |      closed : {'right', 'left'}
 |          Which side of bin interval is closed. The default is 'left'
 |          for all frequency offsets except for 'M', 'A', 'Q', 'BM',
 |          'BA', 'BQ', and 'W' which all have a default of 'right'.
 |      label : {'right', 'left'}
 |          Which bin edge label to label bucket with. The default is 'left'
 |          for all frequency offsets except for 'M', 'A', 'Q', 'BM',
 |          'BA', 'BQ', and 'W' which all have a default of 'right'.
 |      convention : {'start', 'end', 's', 'e'}
 |          For PeriodIndex only, controls whether to use the start or end of
 |          `rule`
 |      kind: {'timestamp', 'period'}, optional
 |          Pass 'timestamp' to convert the resulting index to a
 |          ``DateTimeIndex`` or 'period' to convert it to a ``PeriodIndex``.
 |          By default the input representation is retained.
 |      loffset : timedelta
 |          Adjust the resampled time labels
 |      base : int, default 0
 |          For frequencies that evenly subdivide 1 day, the "origin" of the
 |          aggregated intervals. For example, for '5min' frequency, base could
 |          range from 0 through 4. Defaults to 0
 |      on : string, optional
 |          For a DataFrame, column to use instead of index for resampling.
 |          Column must be datetime-like.
 |      
 |          .. versionadded:: 0.19.0
 |      
 |      level : string or int, optional
 |          For a MultiIndex, level (name or number) to use for
 |          resampling.  Level must be datetime-like.
 |      
 |          .. versionadded:: 0.19.0
 |      
 |      Returns
 |      -------
 |      Resampler object
 |      
 |      Notes
 |      -----
 |      See the `user guide
 |      <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#resampling>`_
 |      for more.
 |      
 |      To learn more about the offset strings, please see `this link
 |      <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__.
 |      
 |      Examples
 |      --------
 |      
 |      Start by creating a series with 9 one minute timestamps.
 |      
 |      >>> index = pd.date_range('1/1/2000', periods=9, freq='T')
 |      >>> series = pd.Series(range(9), index=index)
 |      >>> series
 |      2000-01-01 00:00:00    0
 |      2000-01-01 00:01:00    1
 |      2000-01-01 00:02:00    2
 |      2000-01-01 00:03:00    3
 |      2000-01-01 00:04:00    4
 |      2000-01-01 00:05:00    5
 |      2000-01-01 00:06:00    6
 |      2000-01-01 00:07:00    7
 |      2000-01-01 00:08:00    8
 |      Freq: T, dtype: int64
 |      
 |      Downsample the series into 3 minute bins and sum the values
 |      of the timestamps falling into a bin.
 |      
 |      >>> series.resample('3T').sum()
 |      2000-01-01 00:00:00     3
 |      2000-01-01 00:03:00    12
 |      2000-01-01 00:06:00    21
 |      Freq: 3T, dtype: int64
 |      
 |      Downsample the series into 3 minute bins as above, but label each
 |      bin using the right edge instead of the left. Please note that the
 |      value in the bucket used as the label is not included in the bucket,
 |      which it labels. For example, in the original series the
 |      bucket ``2000-01-01 00:03:00`` contains the value 3, but the summed
 |      value in the resampled bucket with the label ``2000-01-01 00:03:00``
 |      does not include 3 (if it did, the summed value would be 6, not 3).
 |      To include this value close the right side of the bin interval as
 |      illustrated in the example below this one.
 |      
 |      >>> series.resample('3T', label='right').sum()
 |      2000-01-01 00:03:00     3
 |      2000-01-01 00:06:00    12
 |      2000-01-01 00:09:00    21
 |      Freq: 3T, dtype: int64
 |      
 |      Downsample the series into 3 minute bins as above, but close the right
 |      side of the bin interval.
 |      
 |      >>> series.resample('3T', label='right', closed='right').sum()
 |      2000-01-01 00:00:00     0
 |      2000-01-01 00:03:00     6
 |      2000-01-01 00:06:00    15
 |      2000-01-01 00:09:00    15
 |      Freq: 3T, dtype: int64
 |      
 |      Upsample the series into 30 second bins.
 |      
 |      >>> series.resample('30S').asfreq()[0:5] #select first 5 rows
 |      2000-01-01 00:00:00   0.0
 |      2000-01-01 00:00:30   NaN
 |      2000-01-01 00:01:00   1.0
 |      2000-01-01 00:01:30   NaN
 |      2000-01-01 00:02:00   2.0
 |      Freq: 30S, dtype: float64
 |      
 |      Upsample the series into 30 second bins and fill the ``NaN``
 |      values using the ``pad`` method.
 |      
 |      >>> series.resample('30S').pad()[0:5]
 |      2000-01-01 00:00:00    0
 |      2000-01-01 00:00:30    0
 |      2000-01-01 00:01:00    1
 |      2000-01-01 00:01:30    1
 |      2000-01-01 00:02:00    2
 |      Freq: 30S, dtype: int64
 |      
 |      Upsample the series into 30 second bins and fill the
 |      ``NaN`` values using the ``bfill`` method.
 |      
 |      >>> series.resample('30S').bfill()[0:5]
 |      2000-01-01 00:00:00    0
 |      2000-01-01 00:00:30    1
 |      2000-01-01 00:01:00    1
 |      2000-01-01 00:01:30    2
 |      2000-01-01 00:02:00    2
 |      Freq: 30S, dtype: int64
 |      
 |      Pass a custom function via ``apply``
 |      
 |      >>> def custom_resampler(array_like):
 |      ...     return np.sum(array_like)+5
 |      
 |      >>> series.resample('3T').apply(custom_resampler)
 |      2000-01-01 00:00:00     8
 |      2000-01-01 00:03:00    17
 |      2000-01-01 00:06:00    26
 |      Freq: 3T, dtype: int64
 |      
 |      For a Series with a PeriodIndex, the keyword `convention` can be
 |      used to control whether to use the start or end of `rule`.
 |      
 |      >>> s = pd.Series([1, 2], index=pd.period_range('2012-01-01',
 |                                                      freq='A',
 |                                                      periods=2))
 |      >>> s
 |      2012    1
 |      2013    2
 |      Freq: A-DEC, dtype: int64
 |      
 |      Resample by month using 'start' `convention`. Values are assigned to
 |      the first month of the period.
 |      
 |      >>> s.resample('M', convention='start').asfreq().head()
 |      2012-01    1.0
 |      2012-02    NaN
 |      2012-03    NaN
 |      2012-04    NaN
 |      2012-05    NaN
 |      Freq: M, dtype: float64
 |      
 |      Resample by month using 'end' `convention`. Values are assigned to
 |      the last month of the period.
 |      
 |      >>> s.resample('M', convention='end').asfreq()
 |      2012-12    1.0
 |      2013-01    NaN
 |      2013-02    NaN
 |      2013-03    NaN
 |      2013-04    NaN
 |      2013-05    NaN
 |      2013-06    NaN
 |      2013-07    NaN
 |      2013-08    NaN
 |      2013-09    NaN
 |      2013-10    NaN
 |      2013-11    NaN
 |      2013-12    2.0
 |      Freq: M, dtype: float64
 |      
 |      For DataFrame objects, the keyword ``on`` can be used to specify the
 |      column instead of the index for resampling.
 |      
 |      >>> df = pd.DataFrame(data=9*[range(4)], columns=['a', 'b', 'c', 'd'])
 |      >>> df['time'] = pd.date_range('1/1/2000', periods=9, freq='T')
 |      >>> df.resample('3T', on='time').sum()
 |                           a  b  c  d
 |      time
 |      2000-01-01 00:00:00  0  3  6  9
 |      2000-01-01 00:03:00  0  3  6  9
 |      2000-01-01 00:06:00  0  3  6  9
 |      
 |      For a DataFrame with MultiIndex, the keyword ``level`` can be used to
 |      specify on level the resampling needs to take place.
 |      
 |      >>> time = pd.date_range('1/1/2000', periods=5, freq='T')
 |      >>> df2 = pd.DataFrame(data=10*[range(4)],
 |                             columns=['a', 'b', 'c', 'd'],
 |                             index=pd.MultiIndex.from_product([time, [1, 2]])
 |                             )
 |      >>> df2.resample('3T', level=0).sum()
 |                           a  b   c   d
 |      2000-01-01 00:00:00  0  6  12  18
 |      2000-01-01 00:03:00  0  4   8  12
 |      
 |      See also
 |      --------
 |      groupby : Group by mapping, function, label, or list of labels.
 |  
 |  sample(self, n=None, frac=None, replace=False, weights=None, random_state=None, axis=None)
 |      Return a random sample of items from an axis of object.
 |      
 |      You can use `random_state` for reproducibility.
 |      
 |      Parameters
 |      ----------
 |      n : int, optional
 |          Number of items from axis to return. Cannot be used with `frac`.
 |          Default = 1 if `frac` = None.
 |      frac : float, optional
 |          Fraction of axis items to return. Cannot be used with `n`.
 |      replace : boolean, optional
 |          Sample with or without replacement. Default = False.
 |      weights : str or ndarray-like, optional
 |          Default 'None' results in equal probability weighting.
 |          If passed a Series, will align with target object on index. Index
 |          values in weights not found in sampled object will be ignored and
 |          index values in sampled object not in weights will be assigned
 |          weights of zero.
 |          If called on a DataFrame, will accept the name of a column
 |          when axis = 0.
 |          Unless weights are a Series, weights must be same length as axis
 |          being sampled.
 |          If weights do not sum to 1, they will be normalized to sum to 1.
 |          Missing values in the weights column will be treated as zero.
 |          inf and -inf values not allowed.
 |      random_state : int or numpy.random.RandomState, optional
 |          Seed for the random number generator (if int), or numpy RandomState
 |          object.
 |      axis : int or string, optional
 |          Axis to sample. Accepts axis number or name. Default is stat axis
 |          for given data type (0 for Series and DataFrames, 1 for Panels).
 |      
 |      Returns
 |      -------
 |      A new object of same type as caller.
 |      
 |      Examples
 |      --------
 |      Generate an example ``Series`` and ``DataFrame``:
 |      
 |      >>> s = pd.Series(np.random.randn(50))
 |      >>> s.head()
 |      0   -0.038497
 |      1    1.820773
 |      2   -0.972766
 |      3   -1.598270
 |      4   -1.095526
 |      dtype: float64
 |      >>> df = pd.DataFrame(np.random.randn(50, 4), columns=list('ABCD'))
 |      >>> df.head()
 |                A         B         C         D
 |      0  0.016443 -2.318952 -0.566372 -1.028078
 |      1 -1.051921  0.438836  0.658280 -0.175797
 |      2 -1.243569 -0.364626 -0.215065  0.057736
 |      3  1.768216  0.404512 -0.385604 -1.457834
 |      4  1.072446 -1.137172  0.314194 -0.046661
 |      
 |      Next extract a random sample from both of these objects...
 |      
 |      3 random elements from the ``Series``:
 |      
 |      >>> s.sample(n=3)
 |      27   -0.994689
 |      55   -1.049016
 |      67   -0.224565
 |      dtype: float64
 |      
 |      And a random 10% of the ``DataFrame`` with replacement:
 |      
 |      >>> df.sample(frac=0.1, replace=True)
 |                 A         B         C         D
 |      35  1.981780  0.142106  1.817165 -0.290805
 |      49 -1.336199 -0.448634 -0.789640  0.217116
 |      40  0.823173 -0.078816  1.009536  1.015108
 |      15  1.421154 -0.055301 -1.922594 -0.019696
 |      6  -0.148339  0.832938  1.787600 -1.383767
 |      
 |      You can use `random state` for reproducibility:
 |      
 |      >>> df.sample(random_state=1)
 |      A         B         C         D
 |      37 -2.027662  0.103611  0.237496 -0.165867
 |      43 -0.259323 -0.583426  1.516140 -0.479118
 |      12 -1.686325 -0.579510  0.985195 -0.460286
 |      8   1.167946  0.429082  1.215742 -1.636041
 |      9   1.197475 -0.864188  1.554031 -1.505264
 |  
 |  select(self, crit, axis=0)
 |      Return data corresponding to axis labels matching criteria
 |      
 |      .. deprecated:: 0.21.0
 |          Use df.loc[df.index.map(crit)] to select via labels
 |      
 |      Parameters
 |      ----------
 |      crit : function
 |          To be called on each index (label). Should return True or False
 |      axis : int
 |      
 |      Returns
 |      -------
 |      selection : type of caller
 |  
 |  set_axis(self, labels, axis=0, inplace=None)
 |      Assign desired index to given axis.
 |      
 |      Indexes for column or row labels can be changed by assigning
 |      a list-like or Index.
 |      
 |      .. versionchanged:: 0.21.0
 |      
 |         The signature is now `labels` and `axis`, consistent with
 |         the rest of pandas API. Previously, the `axis` and `labels`
 |         arguments were respectively the first and second positional
 |         arguments.
 |      
 |      Parameters
 |      ----------
 |      labels : list-like, Index
 |          The values for the new index.
 |      
 |      axis : {0 or 'index', 1 or 'columns'}, default 0
 |          The axis to update. The value 0 identifies the rows, and 1
 |          identifies the columns.
 |      
 |      inplace : boolean, default None
 |          Whether to return a new %(klass)s instance.
 |      
 |          .. warning::
 |      
 |             ``inplace=None`` currently falls back to to True, but in a
 |             future version, will default to False. Use inplace=True
 |             explicitly rather than relying on the default.
 |      
 |      Returns
 |      -------
 |      renamed : %(klass)s or None
 |          An object of same type as caller if inplace=False, None otherwise.
 |      
 |      See Also
 |      --------
 |      pandas.DataFrame.rename_axis : Alter the name of the index or columns.
 |      
 |      Examples
 |      --------
 |      **Series**
 |      
 |      >>> s = pd.Series([1, 2, 3])
 |      >>> s
 |      0    1
 |      1    2
 |      2    3
 |      dtype: int64
 |      
 |      >>> s.set_axis(['a', 'b', 'c'], axis=0, inplace=False)
 |      a    1
 |      b    2
 |      c    3
 |      dtype: int64
 |      
 |      The original object is not modified.
 |      
 |      >>> s
 |      0    1
 |      1    2
 |      2    3
 |      dtype: int64
 |      
 |      **DataFrame**
 |      
 |      >>> df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
 |      
 |      Change the row labels.
 |      
 |      >>> df.set_axis(['a', 'b', 'c'], axis='index', inplace=False)
 |         A  B
 |      a  1  4
 |      b  2  5
 |      c  3  6
 |      
 |      Change the column labels.
 |      
 |      >>> df.set_axis(['I', 'II'], axis='columns', inplace=False)
 |         I  II
 |      0  1   4
 |      1  2   5
 |      2  3   6
 |      
 |      Now, update the labels inplace.
 |      
 |      >>> df.set_axis(['i', 'ii'], axis='columns', inplace=True)
 |      >>> df
 |         i  ii
 |      0  1   4
 |      1  2   5
 |      2  3   6
 |  
 |  slice_shift(self, periods=1, axis=0)
 |      Equivalent to `shift` without copying data. The shifted data will
 |      not include the dropped periods and the shifted axis will be smaller
 |      than the original.
 |      
 |      Parameters
 |      ----------
 |      periods : int
 |          Number of periods to move, can be positive or negative
 |      
 |      Notes
 |      -----
 |      While the `slice_shift` is faster than `shift`, you may pay for it
 |      later during alignment.
 |      
 |      Returns
 |      -------
 |      shifted : same type as caller
 |  
 |  squeeze(self, axis=None)
 |      Squeeze length 1 dimensions.
 |      
 |      Parameters
 |      ----------
 |      axis : None, integer or string axis name, optional
 |          The axis to squeeze if 1-sized.
 |      
 |          .. versionadded:: 0.20.0
 |      
 |      Returns
 |      -------
 |      scalar if 1-sized, else original object
 |  
 |  swapaxes(self, axis1, axis2, copy=True)
 |      Interchange axes and swap values axes appropriately
 |      
 |      Returns
 |      -------
 |      y : same as input
 |  
 |  tail(self, n=5)
 |      Return the last `n` rows.
 |      
 |      This function returns last `n` rows from the object based on
 |      position. It is useful for quickly verifying data, for example,
 |      after sorting or appending rows.
 |      
 |      Parameters
 |      ----------
 |      n : int, default 5
 |          Number of rows to select.
 |      
 |      Returns
 |      -------
 |      type of caller
 |          The last `n` rows of the caller object.
 |      
 |      See Also
 |      --------
 |      pandas.DataFrame.head : The first `n` rows of the caller object.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'animal':['alligator', 'bee', 'falcon', 'lion',
 |      ...                    'monkey', 'parrot', 'shark', 'whale', 'zebra']})
 |      >>> df
 |            animal
 |      0  alligator
 |      1        bee
 |      2     falcon
 |      3       lion
 |      4     monkey
 |      5     parrot
 |      6      shark
 |      7      whale
 |      8      zebra
 |      
 |      Viewing the last 5 lines
 |      
 |      >>> df.tail()
 |         animal
 |      4  monkey
 |      5  parrot
 |      6   shark
 |      7   whale
 |      8   zebra
 |      
 |      Viewing the last `n` lines (three in this case)
 |      
 |      >>> df.tail(3)
 |        animal
 |      6  shark
 |      7  whale
 |      8  zebra
 |  
 |  take(self, indices, axis=0, convert=None, is_copy=True, **kwargs)
 |      Return the elements in the given *positional* indices along an axis.
 |      
 |      This means that we are not indexing according to actual values in
 |      the index attribute of the object. We are indexing according to the
 |      actual position of the element in the object.
 |      
 |      Parameters
 |      ----------
 |      indices : array-like
 |          An array of ints indicating which positions to take.
 |      axis : {0 or 'index', 1 or 'columns', None}, default 0
 |          The axis on which to select elements. ``0`` means that we are
 |          selecting rows, ``1`` means that we are selecting columns.
 |      convert : bool, default True
 |          Whether to convert negative indices into positive ones.
 |          For example, ``-1`` would map to the ``len(axis) - 1``.
 |          The conversions are similar to the behavior of indexing a
 |          regular Python list.
 |      
 |          .. deprecated:: 0.21.0
 |             In the future, negative indices will always be converted.
 |      
 |      is_copy : bool, default True
 |          Whether to return a copy of the original object or not.
 |      **kwargs
 |          For compatibility with :meth:`numpy.take`. Has no effect on the
 |          output.
 |      
 |      Returns
 |      -------
 |      taken : type of caller
 |          An array-like containing the elements taken from the object.
 |      
 |      See Also
 |      --------
 |      DataFrame.loc : Select a subset of a DataFrame by labels.
 |      DataFrame.iloc : Select a subset of a DataFrame by positions.
 |      numpy.take : Take elements from an array along an axis.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame([('falcon', 'bird',    389.0),
 |      ...                    ('parrot', 'bird',     24.0),
 |      ...                    ('lion',   'mammal',   80.5),
 |      ...                    ('monkey', 'mammal', np.nan)],
 |      ...                    columns=['name', 'class', 'max_speed'],
 |      ...                    index=[0, 2, 3, 1])
 |      >>> df
 |           name   class  max_speed
 |      0  falcon    bird      389.0
 |      2  parrot    bird       24.0
 |      3    lion  mammal       80.5
 |      1  monkey  mammal        NaN
 |      
 |      Take elements at positions 0 and 3 along the axis 0 (default).
 |      
 |      Note how the actual indices selected (0 and 1) do not correspond to
 |      our selected indices 0 and 3. That's because we are selecting the 0th
 |      and 3rd rows, not rows whose indices equal 0 and 3.
 |      
 |      >>> df.take([0, 3])
 |           name   class  max_speed
 |      0  falcon    bird      389.0
 |      1  monkey  mammal        NaN
 |      
 |      Take elements at indices 1 and 2 along the axis 1 (column selection).
 |      
 |      >>> df.take([1, 2], axis=1)
 |          class  max_speed
 |      0    bird      389.0
 |      2    bird       24.0
 |      3  mammal       80.5
 |      1  mammal        NaN
 |      
 |      We may take elements using negative integers for positive indices,
 |      starting from the end of the object, just like with Python lists.
 |      
 |      >>> df.take([-1, -2])
 |           name   class  max_speed
 |      1  monkey  mammal        NaN
 |      3    lion  mammal       80.5
 |  
 |  to_clipboard(self, excel=True, sep=None, **kwargs)
 |      Copy object to the system clipboard.
 |      
 |      Write a text representation of object to the system clipboard.
 |      This can be pasted into Excel, for example.
 |      
 |      Parameters
 |      ----------
 |      excel : bool, default True
 |          - True, use the provided separator, writing in a csv format for
 |            allowing easy pasting into excel.
 |          - False, write a string representation of the object to the
 |            clipboard.
 |      
 |      sep : str, default ``'\t'``
 |          Field delimiter.
 |      **kwargs
 |          These parameters will be passed to DataFrame.to_csv.
 |      
 |      See Also
 |      --------
 |      DataFrame.to_csv : Write a DataFrame to a comma-separated values
 |          (csv) file.
 |      read_clipboard : Read text from clipboard and pass to read_table.
 |      
 |      Notes
 |      -----
 |      Requirements for your platform.
 |      
 |        - Linux : `xclip`, or `xsel` (with `gtk` or `PyQt4` modules)
 |        - Windows : none
 |        - OS X : none
 |      
 |      Examples
 |      --------
 |      Copy the contents of a DataFrame to the clipboard.
 |      
 |      >>> df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=['A', 'B', 'C'])
 |      >>> df.to_clipboard(sep=',')
 |      ... # Wrote the following to the system clipboard:
 |      ... # ,A,B,C
 |      ... # 0,1,2,3
 |      ... # 1,4,5,6
 |      
 |      We can omit the the index by passing the keyword `index` and setting
 |      it to false.
 |      
 |      >>> df.to_clipboard(sep=',', index=False)
 |      ... # Wrote the following to the system clipboard:
 |      ... # A,B,C
 |      ... # 1,2,3
 |      ... # 4,5,6
 |  
 |  to_dense(self)
 |      Return dense representation of NDFrame (as opposed to sparse)
 |  
 |  to_hdf(self, path_or_buf, key, **kwargs)
 |      Write the contained data to an HDF5 file using HDFStore.
 |      
 |      Hierarchical Data Format (HDF) is self-describing, allowing an
 |      application to interpret the structure and contents of a file with
 |      no outside information. One HDF file can hold a mix of related objects
 |      which can be accessed as a group or as individual objects.
 |      
 |      In order to add another DataFrame or Series to an existing HDF file
 |      please use append mode and a different a key.
 |      
 |      For more information see the :ref:`user guide <io.hdf5>`.
 |      
 |      Parameters
 |      ----------
 |      path_or_buf : str or pandas.HDFStore
 |          File path or HDFStore object.
 |      key : str
 |          Identifier for the group in the store.
 |      mode : {'a', 'w', 'r+'}, default 'a'
 |          Mode to open file:
 |      
 |          - 'w': write, a new file is created (an existing file with
 |            the same name would be deleted).
 |          - 'a': append, an existing file is opened for reading and
 |            writing, and if the file does not exist it is created.
 |          - 'r+': similar to 'a', but the file must already exist.
 |      format : {'fixed', 'table'}, default 'fixed'
 |          Possible values:
 |      
 |          - 'fixed': Fixed format. Fast writing/reading. Not-appendable,
 |            nor searchable.
 |          - 'table': Table format. Write as a PyTables Table structure
 |            which may perform worse but allow more flexible operations
 |            like searching / selecting subsets of the data.
 |      append : bool, default False
 |          For Table formats, append the input data to the existing.
 |      data_columns :  list of columns or True, optional
 |          List of columns to create as indexed data columns for on-disk
 |          queries, or True to use all columns. By default only the axes
 |          of the object are indexed. See :ref:`io.hdf5-query-data-columns`.
 |          Applicable only to format='table'.
 |      complevel : {0-9}, optional
 |          Specifies a compression level for data.
 |          A value of 0 disables compression.
 |      complib : {'zlib', 'lzo', 'bzip2', 'blosc'}, default 'zlib'
 |          Specifies the compression library to be used.
 |          As of v0.20.2 these additional compressors for Blosc are supported
 |          (default if no compressor specified: 'blosc:blosclz'):
 |          {'blosc:blosclz', 'blosc:lz4', 'blosc:lz4hc', 'blosc:snappy',
 |          'blosc:zlib', 'blosc:zstd'}.
 |          Specifying a compression library which is not available issues
 |          a ValueError.
 |      fletcher32 : bool, default False
 |          If applying compression use the fletcher32 checksum.
 |      dropna : bool, default False
 |          If true, ALL nan rows will not be written to store.
 |      errors : str, default 'strict'
 |          Specifies how encoding and decoding errors are to be handled.
 |          See the errors argument for :func:`open` for a full list
 |          of options.
 |      
 |      See Also
 |      --------
 |      DataFrame.read_hdf : Read from HDF file.
 |      DataFrame.to_parquet : Write a DataFrame to the binary parquet format.
 |      DataFrame.to_sql : Write to a sql table.
 |      DataFrame.to_feather : Write out feather-format for DataFrames.
 |      DataFrame.to_csv : Write out to a csv file.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]},
 |      ...                   index=['a', 'b', 'c'])
 |      >>> df.to_hdf('data.h5', key='df', mode='w')
 |      
 |      We can add another object to the same file:
 |      
 |      >>> s = pd.Series([1, 2, 3, 4])
 |      >>> s.to_hdf('data.h5', key='s')
 |      
 |      Reading from HDF file:
 |      
 |      >>> pd.read_hdf('data.h5', 'df')
 |      A  B
 |      a  1  4
 |      b  2  5
 |      c  3  6
 |      >>> pd.read_hdf('data.h5', 's')
 |      0    1
 |      1    2
 |      2    3
 |      3    4
 |      dtype: int64
 |      
 |      Deleting file with data:
 |      
 |      >>> import os
 |      >>> os.remove('data.h5')
 |  
 |  to_json(self, path_or_buf=None, orient=None, date_format=None, double_precision=10, force_ascii=True, date_unit='ms', default_handler=None, lines=False, compression=None, index=True)
 |      Convert the object to a JSON string.
 |      
 |      Note NaN's and None will be converted to null and datetime objects
 |      will be converted to UNIX timestamps.
 |      
 |      Parameters
 |      ----------
 |      path_or_buf : string or file handle, optional
 |          File path or object. If not specified, the result is returned as
 |          a string.
 |      orient : string
 |          Indication of expected JSON string format.
 |      
 |          * Series
 |      
 |            - default is 'index'
 |            - allowed values are: {'split','records','index'}
 |      
 |          * DataFrame
 |      
 |            - default is 'columns'
 |            - allowed values are:
 |              {'split','records','index','columns','values'}
 |      
 |          * The format of the JSON string
 |      
 |            - 'split' : dict like {'index' -> [index],
 |              'columns' -> [columns], 'data' -> [values]}
 |            - 'records' : list like
 |              [{column -> value}, ... , {column -> value}]
 |            - 'index' : dict like {index -> {column -> value}}
 |            - 'columns' : dict like {column -> {index -> value}}
 |            - 'values' : just the values array
 |            - 'table' : dict like {'schema': {schema}, 'data': {data}}
 |              describing the data, and the data component is
 |              like ``orient='records'``.
 |      
 |              .. versionchanged:: 0.20.0
 |      
 |      date_format : {None, 'epoch', 'iso'}
 |          Type of date conversion. 'epoch' = epoch milliseconds,
 |          'iso' = ISO8601. The default depends on the `orient`. For
 |          ``orient='table'``, the default is 'iso'. For all other orients,
 |          the default is 'epoch'.
 |      double_precision : int, default 10
 |          The number of decimal places to use when encoding
 |          floating point values.
 |      force_ascii : boolean, default True
 |          Force encoded string to be ASCII.
 |      date_unit : string, default 'ms' (milliseconds)
 |          The time unit to encode to, governs timestamp and ISO8601
 |          precision.  One of 's', 'ms', 'us', 'ns' for second, millisecond,
 |          microsecond, and nanosecond respectively.
 |      default_handler : callable, default None
 |          Handler to call if object cannot otherwise be converted to a
 |          suitable format for JSON. Should receive a single argument which is
 |          the object to convert and return a serialisable object.
 |      lines : boolean, default False
 |          If 'orient' is 'records' write out line delimited json format. Will
 |          throw ValueError if incorrect 'orient' since others are not list
 |          like.
 |      
 |          .. versionadded:: 0.19.0
 |      
 |      compression : {None, 'gzip', 'bz2', 'zip', 'xz'}
 |          A string representing the compression to use in the output file,
 |          only used when the first argument is a filename.
 |      
 |          .. versionadded:: 0.21.0
 |      
 |      index : boolean, default True
 |          Whether to include the index values in the JSON string. Not
 |          including the index (``index=False``) is only supported when
 |          orient is 'split' or 'table'.
 |      
 |          .. versionadded:: 0.23.0
 |      
 |      See Also
 |      --------
 |      pandas.read_json
 |      
 |      Examples
 |      --------
 |      
 |      >>> df = pd.DataFrame([['a', 'b'], ['c', 'd']],
 |      ...                   index=['row 1', 'row 2'],
 |      ...                   columns=['col 1', 'col 2'])
 |      >>> df.to_json(orient='split')
 |      '{"columns":["col 1","col 2"],
 |        "index":["row 1","row 2"],
 |        "data":[["a","b"],["c","d"]]}'
 |      
 |      Encoding/decoding a Dataframe using ``'records'`` formatted JSON.
 |      Note that index labels are not preserved with this encoding.
 |      
 |      >>> df.to_json(orient='records')
 |      '[{"col 1":"a","col 2":"b"},{"col 1":"c","col 2":"d"}]'
 |      
 |      Encoding/decoding a Dataframe using ``'index'`` formatted JSON:
 |      
 |      >>> df.to_json(orient='index')
 |      '{"row 1":{"col 1":"a","col 2":"b"},"row 2":{"col 1":"c","col 2":"d"}}'
 |      
 |      Encoding/decoding a Dataframe using ``'columns'`` formatted JSON:
 |      
 |      >>> df.to_json(orient='columns')
 |      '{"col 1":{"row 1":"a","row 2":"c"},"col 2":{"row 1":"b","row 2":"d"}}'
 |      
 |      Encoding/decoding a Dataframe using ``'values'`` formatted JSON:
 |      
 |      >>> df.to_json(orient='values')
 |      '[["a","b"],["c","d"]]'
 |      
 |      Encoding with Table Schema
 |      
 |      >>> df.to_json(orient='table')
 |      '{"schema": {"fields": [{"name": "index", "type": "string"},
 |                              {"name": "col 1", "type": "string"},
 |                              {"name": "col 2", "type": "string"}],
 |                   "primaryKey": "index",
 |                   "pandas_version": "0.20.0"},
 |        "data": [{"index": "row 1", "col 1": "a", "col 2": "b"},
 |                 {"index": "row 2", "col 1": "c", "col 2": "d"}]}'
 |  
 |  to_latex(self, buf=None, columns=None, col_space=None, header=True, index=True, na_rep='NaN', formatters=None, float_format=None, sparsify=None, index_names=True, bold_rows=False, column_format=None, longtable=None, escape=None, encoding=None, decimal='.', multicolumn=None, multicolumn_format=None, multirow=None)
 |      Render an object to a tabular environment table. You can splice
 |      this into a LaTeX document. Requires \\usepackage{booktabs}.
 |      
 |      .. versionchanged:: 0.20.2
 |         Added to Series
 |      
 |      `to_latex`-specific options:
 |      
 |      bold_rows : boolean, default False
 |          Make the row labels bold in the output
 |      column_format : str, default None
 |          The columns format as specified in `LaTeX table format
 |          <https://en.wikibooks.org/wiki/LaTeX/Tables>`__ e.g 'rcl' for 3
 |          columns
 |      longtable : boolean, default will be read from the pandas config module
 |          Default: False.
 |          Use a longtable environment instead of tabular. Requires adding
 |          a \\usepackage{longtable} to your LaTeX preamble.
 |      escape : boolean, default will be read from the pandas config module
 |          Default: True.
 |          When set to False prevents from escaping latex special
 |          characters in column names.
 |      encoding : str, default None
 |          A string representing the encoding to use in the output file,
 |          defaults to 'ascii' on Python 2 and 'utf-8' on Python 3.
 |      decimal : string, default '.'
 |          Character recognized as decimal separator, e.g. ',' in Europe.
 |      
 |          .. versionadded:: 0.18.0
 |      
 |      multicolumn : boolean, default True
 |          Use \multicolumn to enhance MultiIndex columns.
 |          The default will be read from the config module.
 |      
 |          .. versionadded:: 0.20.0
 |      
 |      multicolumn_format : str, default 'l'
 |          The alignment for multicolumns, similar to `column_format`
 |          The default will be read from the config module.
 |      
 |          .. versionadded:: 0.20.0
 |      
 |      multirow : boolean, default False
 |          Use \multirow to enhance MultiIndex rows.
 |          Requires adding a \\usepackage{multirow} to your LaTeX preamble.
 |          Will print centered labels (instead of top-aligned)
 |          across the contained rows, separating groups via clines.
 |          The default will be read from the pandas config module.
 |      
 |          .. versionadded:: 0.20.0
 |  
 |  to_msgpack(self, path_or_buf=None, encoding='utf-8', **kwargs)
 |      msgpack (serialize) object to input file path
 |      
 |      THIS IS AN EXPERIMENTAL LIBRARY and the storage format
 |      may not be stable until a future release.
 |      
 |      Parameters
 |      ----------
 |      path : string File path, buffer-like, or None
 |          if None, return generated string
 |      append : boolean whether to append to an existing msgpack
 |          (default is False)
 |      compress : type of compressor (zlib or blosc), default to None (no
 |          compression)
 |  
 |  to_pickle(self, path, compression='infer', protocol=4)
 |      Pickle (serialize) object to file.
 |      
 |      Parameters
 |      ----------
 |      path : str
 |          File path where the pickled object will be stored.
 |      compression : {'infer', 'gzip', 'bz2', 'zip', 'xz', None},         default 'infer'
 |          A string representing the compression to use in the output file. By
 |          default, infers from the file extension in specified path.
 |      
 |          .. versionadded:: 0.20.0
 |      protocol : int
 |          Int which indicates which protocol should be used by the pickler,
 |          default HIGHEST_PROTOCOL (see [1]_ paragraph 12.1.2). The possible
 |          values for this parameter depend on the version of Python. For
 |          Python 2.x, possible values are 0, 1, 2. For Python>=3.0, 3 is a
 |          valid value. For Python >= 3.4, 4 is a valid value. A negative
 |          value for the protocol parameter is equivalent to setting its value
 |          to HIGHEST_PROTOCOL.
 |      
 |          .. [1] https://docs.python.org/3/library/pickle.html
 |          .. versionadded:: 0.21.0
 |      
 |      See Also
 |      --------
 |      read_pickle : Load pickled pandas object (or any object) from file.
 |      DataFrame.to_hdf : Write DataFrame to an HDF5 file.
 |      DataFrame.to_sql : Write DataFrame to a SQL database.
 |      DataFrame.to_parquet : Write a DataFrame to the binary parquet format.
 |      
 |      Examples
 |      --------
 |      >>> original_df = pd.DataFrame({"foo": range(5), "bar": range(5, 10)})
 |      >>> original_df
 |         foo  bar
 |      0    0    5
 |      1    1    6
 |      2    2    7
 |      3    3    8
 |      4    4    9
 |      >>> original_df.to_pickle("./dummy.pkl")
 |      
 |      >>> unpickled_df = pd.read_pickle("./dummy.pkl")
 |      >>> unpickled_df
 |         foo  bar
 |      0    0    5
 |      1    1    6
 |      2    2    7
 |      3    3    8
 |      4    4    9
 |      
 |      >>> import os
 |      >>> os.remove("./dummy.pkl")
 |  
 |  to_sql(self, name, con, schema=None, if_exists='fail', index=True, index_label=None, chunksize=None, dtype=None)
 |      Write records stored in a DataFrame to a SQL database.
 |      
 |      Databases supported by SQLAlchemy [1]_ are supported. Tables can be
 |      newly created, appended to, or overwritten.
 |      
 |      Parameters
 |      ----------
 |      name : string
 |          Name of SQL table.
 |      con : sqlalchemy.engine.Engine or sqlite3.Connection
 |          Using SQLAlchemy makes it possible to use any DB supported by that
 |          library. Legacy support is provided for sqlite3.Connection objects.
 |      schema : string, optional
 |          Specify the schema (if database flavor supports this). If None, use
 |          default schema.
 |      if_exists : {'fail', 'replace', 'append'}, default 'fail'
 |          How to behave if the table already exists.
 |      
 |          * fail: Raise a ValueError.
 |          * replace: Drop the table before inserting new values.
 |          * append: Insert new values to the existing table.
 |      
 |      index : boolean, default True
 |          Write DataFrame index as a column. Uses `index_label` as the column
 |          name in the table.
 |      index_label : string or sequence, default None
 |          Column label for index column(s). If None is given (default) and
 |          `index` is True, then the index names are used.
 |          A sequence should be given if the DataFrame uses MultiIndex.
 |      chunksize : int, optional
 |          Rows will be written in batches of this size at a time. By default,
 |          all rows will be written at once.
 |      dtype : dict, optional
 |          Specifying the datatype for columns. The keys should be the column
 |          names and the values should be the SQLAlchemy types or strings for
 |          the sqlite3 legacy mode.
 |      
 |      Raises
 |      ------
 |      ValueError
 |          When the table already exists and `if_exists` is 'fail' (the
 |          default).
 |      
 |      See Also
 |      --------
 |      pandas.read_sql : read a DataFrame from a table
 |      
 |      References
 |      ----------
 |      .. [1] http://docs.sqlalchemy.org
 |      .. [2] https://www.python.org/dev/peps/pep-0249/
 |      
 |      Examples
 |      --------
 |      
 |      Create an in-memory SQLite database.
 |      
 |      >>> from sqlalchemy import create_engine
 |      >>> engine = create_engine('sqlite://', echo=False)
 |      
 |      Create a table from scratch with 3 rows.
 |      
 |      >>> df = pd.DataFrame({'name' : ['User 1', 'User 2', 'User 3']})
 |      >>> df
 |           name
 |      0  User 1
 |      1  User 2
 |      2  User 3
 |      
 |      >>> df.to_sql('users', con=engine)
 |      >>> engine.execute("SELECT * FROM users").fetchall()
 |      [(0, 'User 1'), (1, 'User 2'), (2, 'User 3')]
 |      
 |      >>> df1 = pd.DataFrame({'name' : ['User 4', 'User 5']})
 |      >>> df1.to_sql('users', con=engine, if_exists='append')
 |      >>> engine.execute("SELECT * FROM users").fetchall()
 |      [(0, 'User 1'), (1, 'User 2'), (2, 'User 3'),
 |       (0, 'User 4'), (1, 'User 5')]
 |      
 |      Overwrite the table with just ``df1``.
 |      
 |      >>> df1.to_sql('users', con=engine, if_exists='replace',
 |      ...            index_label='id')
 |      >>> engine.execute("SELECT * FROM users").fetchall()
 |      [(0, 'User 4'), (1, 'User 5')]
 |      
 |      Specify the dtype (especially useful for integers with missing values).
 |      Notice that while pandas is forced to store the data as floating point,
 |      the database supports nullable integers. When fetching the data with
 |      Python, we get back integer scalars.
 |      
 |      >>> df = pd.DataFrame({"A": [1, None, 2]})
 |      >>> df
 |           A
 |      0  1.0
 |      1  NaN
 |      2  2.0
 |      
 |      >>> from sqlalchemy.types import Integer
 |      >>> df.to_sql('integers', con=engine, index=False,
 |      ...           dtype={"A": Integer()})
 |      
 |      >>> engine.execute("SELECT * FROM integers").fetchall()
 |      [(1,), (None,), (2,)]
 |  
 |  to_xarray(self)
 |      Return an xarray object from the pandas object.
 |      
 |      Returns
 |      -------
 |      a DataArray for a Series
 |      a Dataset for a DataFrame
 |      a DataArray for higher dims
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'A' : [1, 1, 2],
 |                             'B' : ['foo', 'bar', 'foo'],
 |                             'C' : np.arange(4.,7)})
 |      >>> df
 |         A    B    C
 |      0  1  foo  4.0
 |      1  1  bar  5.0
 |      2  2  foo  6.0
 |      
 |      >>> df.to_xarray()
 |      <xarray.Dataset>
 |      Dimensions:  (index: 3)
 |      Coordinates:
 |        * index    (index) int64 0 1 2
 |      Data variables:
 |          A        (index) int64 1 1 2
 |          B        (index) object 'foo' 'bar' 'foo'
 |          C        (index) float64 4.0 5.0 6.0
 |      
 |      >>> df = pd.DataFrame({'A' : [1, 1, 2],
 |                             'B' : ['foo', 'bar', 'foo'],
 |                             'C' : np.arange(4.,7)}
 |                           ).set_index(['B','A'])
 |      >>> df
 |               C
 |      B   A
 |      foo 1  4.0
 |      bar 1  5.0
 |      foo 2  6.0
 |      
 |      >>> df.to_xarray()
 |      <xarray.Dataset>
 |      Dimensions:  (A: 2, B: 2)
 |      Coordinates:
 |        * B        (B) object 'bar' 'foo'
 |        * A        (A) int64 1 2
 |      Data variables:
 |          C        (B, A) float64 5.0 nan 4.0 6.0
 |      
 |      >>> p = pd.Panel(np.arange(24).reshape(4,3,2),
 |                       items=list('ABCD'),
 |                       major_axis=pd.date_range('20130101', periods=3),
 |                       minor_axis=['first', 'second'])
 |      >>> p
 |      <class 'pandas.core.panel.Panel'>
 |      Dimensions: 4 (items) x 3 (major_axis) x 2 (minor_axis)
 |      Items axis: A to D
 |      Major_axis axis: 2013-01-01 00:00:00 to 2013-01-03 00:00:00
 |      Minor_axis axis: first to second
 |      
 |      >>> p.to_xarray()
 |      <xarray.DataArray (items: 4, major_axis: 3, minor_axis: 2)>
 |      array([[[ 0,  1],
 |              [ 2,  3],
 |              [ 4,  5]],
 |             [[ 6,  7],
 |              [ 8,  9],
 |              [10, 11]],
 |             [[12, 13],
 |              [14, 15],
 |              [16, 17]],
 |             [[18, 19],
 |              [20, 21],
 |              [22, 23]]])
 |      Coordinates:
 |        * items       (items) object 'A' 'B' 'C' 'D'
 |        * major_axis  (major_axis) datetime64[ns] 2013-01-01 2013-01-02 2013-01-03  # noqa
 |        * minor_axis  (minor_axis) object 'first' 'second'
 |      
 |      Notes
 |      -----
 |      See the `xarray docs <http://xarray.pydata.org/en/stable/>`__
 |  
 |  truncate(self, before=None, after=None, axis=None, copy=True)
 |      Truncate a Series or DataFrame before and after some index value.
 |      
 |      This is a useful shorthand for boolean indexing based on index
 |      values above or below certain thresholds.
 |      
 |      Parameters
 |      ----------
 |      before : date, string, int
 |          Truncate all rows before this index value.
 |      after : date, string, int
 |          Truncate all rows after this index value.
 |      axis : {0 or 'index', 1 or 'columns'}, optional
 |          Axis to truncate. Truncates the index (rows) by default.
 |      copy : boolean, default is True,
 |          Return a copy of the truncated section.
 |      
 |      Returns
 |      -------
 |      type of caller
 |          The truncated Series or DataFrame.
 |      
 |      See Also
 |      --------
 |      DataFrame.loc : Select a subset of a DataFrame by label.
 |      DataFrame.iloc : Select a subset of a DataFrame by position.
 |      
 |      Notes
 |      -----
 |      If the index being truncated contains only datetime values,
 |      `before` and `after` may be specified as strings instead of
 |      Timestamps.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'A': ['a', 'b', 'c', 'd', 'e'],
 |      ...                    'B': ['f', 'g', 'h', 'i', 'j'],
 |      ...                    'C': ['k', 'l', 'm', 'n', 'o']},
 |      ...                    index=[1, 2, 3, 4, 5])
 |      >>> df
 |         A  B  C
 |      1  a  f  k
 |      2  b  g  l
 |      3  c  h  m
 |      4  d  i  n
 |      5  e  j  o
 |      
 |      >>> df.truncate(before=2, after=4)
 |         A  B  C
 |      2  b  g  l
 |      3  c  h  m
 |      4  d  i  n
 |      
 |      The columns of a DataFrame can be truncated.
 |      
 |      >>> df.truncate(before="A", after="B", axis="columns")
 |         A  B
 |      1  a  f
 |      2  b  g
 |      3  c  h
 |      4  d  i
 |      5  e  j
 |      
 |      For Series, only rows can be truncated.
 |      
 |      >>> df['A'].truncate(before=2, after=4)
 |      2    b
 |      3    c
 |      4    d
 |      Name: A, dtype: object
 |      
 |      The index values in ``truncate`` can be datetimes or string
 |      dates.
 |      
 |      >>> dates = pd.date_range('2016-01-01', '2016-02-01', freq='s')
 |      >>> df = pd.DataFrame(index=dates, data={'A': 1})
 |      >>> df.tail()
 |                           A
 |      2016-01-31 23:59:56  1
 |      2016-01-31 23:59:57  1
 |      2016-01-31 23:59:58  1
 |      2016-01-31 23:59:59  1
 |      2016-02-01 00:00:00  1
 |      
 |      >>> df.truncate(before=pd.Timestamp('2016-01-05'),
 |      ...             after=pd.Timestamp('2016-01-10')).tail()
 |                           A
 |      2016-01-09 23:59:56  1
 |      2016-01-09 23:59:57  1
 |      2016-01-09 23:59:58  1
 |      2016-01-09 23:59:59  1
 |      2016-01-10 00:00:00  1
 |      
 |      Because the index is a DatetimeIndex containing only dates, we can
 |      specify `before` and `after` as strings. They will be coerced to
 |      Timestamps before truncation.
 |      
 |      >>> df.truncate('2016-01-05', '2016-01-10').tail()
 |                           A
 |      2016-01-09 23:59:56  1
 |      2016-01-09 23:59:57  1
 |      2016-01-09 23:59:58  1
 |      2016-01-09 23:59:59  1
 |      2016-01-10 00:00:00  1
 |      
 |      Note that ``truncate`` assumes a 0 value for any unspecified time
 |      component (midnight). This differs from partial string slicing, which
 |      returns any partially matching dates.
 |      
 |      >>> df.loc['2016-01-05':'2016-01-10', :].tail()
 |                           A
 |      2016-01-10 23:59:55  1
 |      2016-01-10 23:59:56  1
 |      2016-01-10 23:59:57  1
 |      2016-01-10 23:59:58  1
 |      2016-01-10 23:59:59  1
 |  
 |  tshift(self, periods=1, freq=None, axis=0)
 |      Shift the time index, using the index's frequency if available.
 |      
 |      Parameters
 |      ----------
 |      periods : int
 |          Number of periods to move, can be positive or negative
 |      freq : DateOffset, timedelta, or time rule string, default None
 |          Increment to use from the tseries module or time rule (e.g. 'EOM')
 |      axis : int or basestring
 |          Corresponds to the axis that contains the Index
 |      
 |      Notes
 |      -----
 |      If freq is not specified then tries to use the freq or inferred_freq
 |      attributes of the index. If neither of those attributes exist, a
 |      ValueError is thrown
 |      
 |      Returns
 |      -------
 |      shifted : NDFrame
 |  
 |  tz_convert(self, tz, axis=0, level=None, copy=True)
 |      Convert tz-aware axis to target time zone.
 |      
 |      Parameters
 |      ----------
 |      tz : string or pytz.timezone object
 |      axis : the axis to convert
 |      level : int, str, default None
 |          If axis ia a MultiIndex, convert a specific level. Otherwise
 |          must be None
 |      copy : boolean, default True
 |          Also make a copy of the underlying data
 |      
 |      Returns
 |      -------
 |      
 |      Raises
 |      ------
 |      TypeError
 |          If the axis is tz-naive.
 |  
 |  tz_localize(self, tz, axis=0, level=None, copy=True, ambiguous='raise')
 |      Localize tz-naive TimeSeries to target time zone.
 |      
 |      Parameters
 |      ----------
 |      tz : string or pytz.timezone object
 |      axis : the axis to localize
 |      level : int, str, default None
 |          If axis ia a MultiIndex, localize a specific level. Otherwise
 |          must be None
 |      copy : boolean, default True
 |          Also make a copy of the underlying data
 |      ambiguous : 'infer', bool-ndarray, 'NaT', default 'raise'
 |          - 'infer' will attempt to infer fall dst-transition hours based on
 |            order
 |          - bool-ndarray where True signifies a DST time, False designates
 |            a non-DST time (note that this flag is only applicable for
 |            ambiguous times)
 |          - 'NaT' will return NaT where there are ambiguous times
 |          - 'raise' will raise an AmbiguousTimeError if there are ambiguous
 |            times
 |      
 |      Returns
 |      -------
 |      
 |      Raises
 |      ------
 |      TypeError
 |          If the TimeSeries is tz-aware and tz is not None.
 |  
 |  where(self, cond, other=nan, inplace=False, axis=None, level=None, errors='raise', try_cast=False, raise_on_error=None)
 |      Return an object of same shape as self and whose corresponding
 |      entries are from self where `cond` is True and otherwise are from
 |      `other`.
 |      
 |      Parameters
 |      ----------
 |      cond : boolean NDFrame, array-like, or callable
 |          Where `cond` is True, keep the original value. Where
 |          False, replace with corresponding value from `other`.
 |          If `cond` is callable, it is computed on the NDFrame and
 |          should return boolean NDFrame or array. The callable must
 |          not change input NDFrame (though pandas doesn't check it).
 |      
 |          .. versionadded:: 0.18.1
 |              A callable can be used as cond.
 |      
 |      other : scalar, NDFrame, or callable
 |          Entries where `cond` is False are replaced with
 |          corresponding value from `other`.
 |          If other is callable, it is computed on the NDFrame and
 |          should return scalar or NDFrame. The callable must not
 |          change input NDFrame (though pandas doesn't check it).
 |      
 |          .. versionadded:: 0.18.1
 |              A callable can be used as other.
 |      
 |      inplace : boolean, default False
 |          Whether to perform the operation in place on the data
 |      axis : alignment axis if needed, default None
 |      level : alignment level if needed, default None
 |      errors : str, {'raise', 'ignore'}, default 'raise'
 |          - ``raise`` : allow exceptions to be raised
 |          - ``ignore`` : suppress exceptions. On error return original object
 |      
 |          Note that currently this parameter won't affect
 |          the results and will always coerce to a suitable dtype.
 |      
 |      try_cast : boolean, default False
 |          try to cast the result back to the input type (if possible),
 |      raise_on_error : boolean, default True
 |          Whether to raise on invalid data types (e.g. trying to where on
 |          strings)
 |      
 |          .. deprecated:: 0.21.0
 |      
 |      Returns
 |      -------
 |      wh : same type as caller
 |      
 |      Notes
 |      -----
 |      The where method is an application of the if-then idiom. For each
 |      element in the calling DataFrame, if ``cond`` is ``True`` the
 |      element is used; otherwise the corresponding element from the DataFrame
 |      ``other`` is used.
 |      
 |      The signature for :func:`DataFrame.where` differs from
 |      :func:`numpy.where`. Roughly ``df1.where(m, df2)`` is equivalent to
 |      ``np.where(m, df1, df2)``.
 |      
 |      For further details and examples see the ``where`` documentation in
 |      :ref:`indexing <indexing.where_mask>`.
 |      
 |      Examples
 |      --------
 |      >>> s = pd.Series(range(5))
 |      >>> s.where(s > 0)
 |      0    NaN
 |      1    1.0
 |      2    2.0
 |      3    3.0
 |      4    4.0
 |      
 |      >>> s.mask(s > 0)
 |      0    0.0
 |      1    NaN
 |      2    NaN
 |      3    NaN
 |      4    NaN
 |      
 |      >>> s.where(s > 1, 10)
 |      0    10.0
 |      1    10.0
 |      2    2.0
 |      3    3.0
 |      4    4.0
 |      
 |      >>> df = pd.DataFrame(np.arange(10).reshape(-1, 2), columns=['A', 'B'])
 |      >>> m = df % 3 == 0
 |      >>> df.where(m, -df)
 |         A  B
 |      0  0 -1
 |      1 -2  3
 |      2 -4 -5
 |      3  6 -7
 |      4 -8  9
 |      >>> df.where(m, -df) == np.where(m, df, -df)
 |            A     B
 |      0  True  True
 |      1  True  True
 |      2  True  True
 |      3  True  True
 |      4  True  True
 |      >>> df.where(m, -df) == df.mask(~m, -df)
 |            A     B
 |      0  True  True
 |      1  True  True
 |      2  True  True
 |      3  True  True
 |      4  True  True
 |      
 |      See Also
 |      --------
 |      :func:`DataFrame.mask`
 |  
 |  xs(self, key, axis=0, level=None, drop_level=True)
 |      Returns a cross-section (row(s) or column(s)) from the
 |      Series/DataFrame. Defaults to cross-section on the rows (axis=0).
 |      
 |      Parameters
 |      ----------
 |      key : object
 |          Some label contained in the index, or partially in a MultiIndex
 |      axis : int, default 0
 |          Axis to retrieve cross-section on
 |      level : object, defaults to first n levels (n=1 or len(key))
 |          In case of a key partially contained in a MultiIndex, indicate
 |          which levels are used. Levels can be referred by label or position.
 |      drop_level : boolean, default True
 |          If False, returns object with same levels as self.
 |      
 |      Examples
 |      --------
 |      >>> df
 |         A  B  C
 |      a  4  5  2
 |      b  4  0  9
 |      c  9  7  3
 |      >>> df.xs('a')
 |      A    4
 |      B    5
 |      C    2
 |      Name: a
 |      >>> df.xs('C', axis=1)
 |      a    2
 |      b    9
 |      c    3
 |      Name: C
 |      
 |      >>> df
 |                          A  B  C  D
 |      first second third
 |      bar   one    1      4  1  8  9
 |            two    1      7  5  5  0
 |      baz   one    1      6  6  8  0
 |            three  2      5  3  5  3
 |      >>> df.xs(('baz', 'three'))
 |             A  B  C  D
 |      third
 |      2      5  3  5  3
 |      >>> df.xs('one', level=1)
 |                   A  B  C  D
 |      first third
 |      bar   1      4  1  8  9
 |      baz   1      6  6  8  0
 |      >>> df.xs(('baz', 2), level=[0, 'third'])
 |              A  B  C  D
 |      second
 |      three   5  3  5  3
 |      
 |      Returns
 |      -------
 |      xs : Series or DataFrame
 |      
 |      Notes
 |      -----
 |      xs is only for getting, not setting values.
 |      
 |      MultiIndex Slicers is a generic way to get/set values on any level or
 |      levels.  It is a superset of xs functionality, see
 |      :ref:`MultiIndex Slicers <advanced.mi_slicers>`
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from pandas.core.generic.NDFrame:
 |  
 |  at
 |      Access a single value for a row/column label pair.
 |      
 |      Similar to ``loc``, in that both provide label-based lookups. Use
 |      ``at`` if you only need to get or set a single value in a DataFrame
 |      or Series.
 |      
 |      See Also
 |      --------
 |      DataFrame.iat : Access a single value for a row/column pair by integer
 |          position
 |      DataFrame.loc : Access a group of rows and columns by label(s)
 |      Series.at : Access a single value using a label
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame([[0, 2, 3], [0, 4, 1], [10, 20, 30]],
 |      ...                   index=[4, 5, 6], columns=['A', 'B', 'C'])
 |      >>> df
 |          A   B   C
 |      4   0   2   3
 |      5   0   4   1
 |      6  10  20  30
 |      
 |      Get value at specified row/column pair
 |      
 |      >>> df.at[4, 'B']
 |      2
 |      
 |      Set value at specified row/column pair
 |      
 |      >>> df.at[4, 'B'] = 10
 |      >>> df.at[4, 'B']
 |      10
 |      
 |      Get value within a Series
 |      
 |      >>> df.loc[5].at['B']
 |      4
 |      
 |      Raises
 |      ------
 |      KeyError
 |          When label does not exist in DataFrame
 |  
 |  blocks
 |      Internal property, property synonym for as_blocks()
 |      
 |      .. deprecated:: 0.21.0
 |  
 |  dtypes
 |      Return the dtypes in the DataFrame.
 |      
 |      This returns a Series with the data type of each column.
 |      The result's index is the original DataFrame's columns. Columns
 |      with mixed types are stored with the ``object`` dtype. See
 |      :ref:`the User Guide <basics.dtypes>` for more.
 |      
 |      Returns
 |      -------
 |      pandas.Series
 |          The data type of each column.
 |      
 |      See Also
 |      --------
 |      pandas.DataFrame.ftypes : dtype and sparsity information.
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame({'float': [1.0],
 |      ...                    'int': [1],
 |      ...                    'datetime': [pd.Timestamp('20180310')],
 |      ...                    'string': ['foo']})
 |      >>> df.dtypes
 |      float              float64
 |      int                  int64
 |      datetime    datetime64[ns]
 |      string              object
 |      dtype: object
 |  
 |  empty
 |      Indicator whether DataFrame is empty.
 |      
 |      True if DataFrame is entirely empty (no items), meaning any of the
 |      axes are of length 0.
 |      
 |      Returns
 |      -------
 |      bool
 |          If DataFrame is empty, return True, if not return False.
 |      
 |      Notes
 |      -----
 |      If DataFrame contains only NaNs, it is still not considered empty. See
 |      the example below.
 |      
 |      Examples
 |      --------
 |      An example of an actual empty DataFrame. Notice the index is empty:
 |      
 |      >>> df_empty = pd.DataFrame({'A' : []})
 |      >>> df_empty
 |      Empty DataFrame
 |      Columns: [A]
 |      Index: []
 |      >>> df_empty.empty
 |      True
 |      
 |      If we only have NaNs in our DataFrame, it is not considered empty! We
 |      will need to drop the NaNs to make the DataFrame empty:
 |      
 |      >>> df = pd.DataFrame({'A' : [np.nan]})
 |      >>> df
 |          A
 |      0 NaN
 |      >>> df.empty
 |      False
 |      >>> df.dropna().empty
 |      True
 |      
 |      See also
 |      --------
 |      pandas.Series.dropna
 |      pandas.DataFrame.dropna
 |  
 |  ftypes
 |      Return the ftypes (indication of sparse/dense and dtype) in DataFrame.
 |      
 |      This returns a Series with the data type of each column.
 |      The result's index is the original DataFrame's columns. Columns
 |      with mixed types are stored with the ``object`` dtype.  See
 |      :ref:`the User Guide <basics.dtypes>` for more.
 |      
 |      Returns
 |      -------
 |      pandas.Series
 |          The data type and indication of sparse/dense of each column.
 |      
 |      See Also
 |      --------
 |      pandas.DataFrame.dtypes: Series with just dtype information.
 |      pandas.SparseDataFrame : Container for sparse tabular data.
 |      
 |      Notes
 |      -----
 |      Sparse data should have the same dtypes as its dense representation.
 |      
 |      Examples
 |      --------
 |      >>> import numpy as np
 |      >>> arr = np.random.RandomState(0).randn(100, 4)
 |      >>> arr[arr < .8] = np.nan
 |      >>> pd.DataFrame(arr).ftypes
 |      0    float64:dense
 |      1    float64:dense
 |      2    float64:dense
 |      3    float64:dense
 |      dtype: object
 |      
 |      >>> pd.SparseDataFrame(arr).ftypes
 |      0    float64:sparse
 |      1    float64:sparse
 |      2    float64:sparse
 |      3    float64:sparse
 |      dtype: object
 |  
 |  iat
 |      Access a single value for a row/column pair by integer position.
 |      
 |      Similar to ``iloc``, in that both provide integer-based lookups. Use
 |      ``iat`` if you only need to get or set a single value in a DataFrame
 |      or Series.
 |      
 |      See Also
 |      --------
 |      DataFrame.at : Access a single value for a row/column label pair
 |      DataFrame.loc : Access a group of rows and columns by label(s)
 |      DataFrame.iloc : Access a group of rows and columns by integer position(s)
 |      
 |      Examples
 |      --------
 |      >>> df = pd.DataFrame([[0, 2, 3], [0, 4, 1], [10, 20, 30]],
 |      ...                   columns=['A', 'B', 'C'])
 |      >>> df
 |          A   B   C
 |      0   0   2   3
 |      1   0   4   1
 |      2  10  20  30
 |      
 |      Get value at specified row/column pair
 |      
 |      >>> df.iat[1, 2]
 |      1
 |      
 |      Set value at specified row/column pair
 |      
 |      >>> df.iat[1, 2] = 10
 |      >>> df.iat[1, 2]
 |      10
 |      
 |      Get value within a series
 |      
 |      >>> df.loc[0].iat[1]
 |      2
 |      
 |      Raises
 |      ------
 |      IndexError
 |          When integer position is out of bounds
 |  
 |  iloc
 |      Purely integer-location based indexing for selection by position.
 |      
 |      ``.iloc[]`` is primarily integer position based (from ``0`` to
 |      ``length-1`` of the axis), but may also be used with a boolean
 |      array.
 |      
 |      Allowed inputs are:
 |      
 |      - An integer, e.g. ``5``.
 |      - A list or array of integers, e.g. ``[4, 3, 0]``.
 |      - A slice object with ints, e.g. ``1:7``.
 |      - A boolean array.
 |      - A ``callable`` function with one argument (the calling Series, DataFrame
 |        or Panel) and that returns valid output for indexing (one of the above)
 |      
 |      ``.iloc`` will raise ``IndexError`` if a requested indexer is
 |      out-of-bounds, except *slice* indexers which allow out-of-bounds
 |      indexing (this conforms with python/numpy *slice* semantics).
 |      
 |      See more at :ref:`Selection by Position <indexing.integer>`
 |  
 |  is_copy
 |  
 |  ix
 |      A primarily label-location based indexer, with integer position
 |      fallback.
 |      
 |      Warning: Starting in 0.20.0, the .ix indexer is deprecated, in
 |      favor of the more strict .iloc and .loc indexers.
 |      
 |      ``.ix[]`` supports mixed integer and label based access. It is
 |      primarily label based, but will fall back to integer positional
 |      access unless the corresponding axis is of integer type.
 |      
 |      ``.ix`` is the most general indexer and will support any of the
 |      inputs in ``.loc`` and ``.iloc``. ``.ix`` also supports floating
 |      point label schemes. ``.ix`` is exceptionally useful when dealing
 |      with mixed positional and label based hierarchical indexes.
 |      
 |      However, when an axis is integer based, ONLY label based access
 |      and not positional access is supported. Thus, in such cases, it's
 |      usually better to be explicit and use ``.iloc`` or ``.loc``.
 |      
 |      See more at :ref:`Advanced Indexing <advanced>`.
 |  
 |  loc
 |      Access a group of rows and columns by label(s) or a boolean array.
 |      
 |      ``.loc[]`` is primarily label based, but may also be used with a
 |      boolean array.
 |      
 |      Allowed inputs are:
 |      
 |      - A single label, e.g. ``5`` or ``'a'``, (note that ``5`` is
 |        interpreted as a *label* of the index, and **never** as an
 |        integer position along the index).
 |      - A list or array of labels, e.g. ``['a', 'b', 'c']``.
 |      - A slice object with labels, e.g. ``'a':'f'``.
 |      
 |        .. warning:: Note that contrary to usual python slices, **both** the
 |            start and the stop are included
 |      
 |      - A boolean array of the same length as the axis being sliced,
 |        e.g. ``[True, False, True]``.
 |      - A ``callable`` function with one argument (the calling Series, DataFrame
 |        or Panel) and that returns valid output for indexing (one of the above)
 |      
 |      See more at :ref:`Selection by Label <indexing.label>`
 |      
 |      See Also
 |      --------
 |      DataFrame.at : Access a single value for a row/column label pair
 |      DataFrame.iloc : Access group of rows and columns by integer position(s)
 |      DataFrame.xs : Returns a cross-section (row(s) or column(s)) from the
 |          Series/DataFrame.
 |      Series.loc : Access group of values using labels
 |      
 |      Examples
 |      --------
 |      **Getting values**
 |      
 |      >>> df = pd.DataFrame([[1, 2], [4, 5], [7, 8]],
 |      ...      index=['cobra', 'viper', 'sidewinder'],
 |      ...      columns=['max_speed', 'shield'])
 |      >>> df
 |                  max_speed  shield
 |      cobra               1       2
 |      viper               4       5
 |      sidewinder          7       8
 |      
 |      Single label. Note this returns the row as a Series.
 |      
 |      >>> df.loc['viper']
 |      max_speed    4
 |      shield       5
 |      Name: viper, dtype: int64
 |      
 |      List of labels. Note using ``[[]]`` returns a DataFrame.
 |      
 |      >>> df.loc[['viper', 'sidewinder']]
 |                  max_speed  shield
 |      viper               4       5
 |      sidewinder          7       8
 |      
 |      Single label for row and column
 |      
 |      >>> df.loc['cobra', 'shield']
 |      2
 |      
 |      Slice with labels for row and single label for column. As mentioned
 |      above, note that both the start and stop of the slice are included.
 |      
 |      >>> df.loc['cobra':'viper', 'max_speed']
 |      cobra    1
 |      viper    4
 |      Name: max_speed, dtype: int64
 |      
 |      Boolean list with the same length as the row axis
 |      
 |      >>> df.loc[[False, False, True]]
 |                  max_speed  shield
 |      sidewinder          7       8
 |      
 |      Conditional that returns a boolean Series
 |      
 |      >>> df.loc[df['shield'] > 6]
 |                  max_speed  shield
 |      sidewinder          7       8
 |      
 |      Conditional that returns a boolean Series with column labels specified
 |      
 |      >>> df.loc[df['shield'] > 6, ['max_speed']]
 |                  max_speed
 |      sidewinder          7
 |      
 |      Callable that returns a boolean Series
 |      
 |      >>> df.loc[lambda df: df['shield'] == 8]
 |                  max_speed  shield
 |      sidewinder          7       8
 |      
 |      **Setting values**
 |      
 |      Set value for all items matching the list of labels
 |      
 |      >>> df.loc[['viper', 'sidewinder'], ['shield']] = 50
 |      >>> df
 |                  max_speed  shield
 |      cobra               1       2
 |      viper               4      50
 |      sidewinder          7      50
 |      
 |      Set value for an entire row
 |      
 |      >>> df.loc['cobra'] = 10
 |      >>> df
 |                  max_speed  shield
 |      cobra              10      10
 |      viper               4      50
 |      sidewinder          7      50
 |      
 |      Set value for an entire column
 |      
 |      >>> df.loc[:, 'max_speed'] = 30
 |      >>> df
 |                  max_speed  shield
 |      cobra              30      10
 |      viper              30      50
 |      sidewinder         30      50
 |      
 |      Set value for rows matching callable condition
 |      
 |      >>> df.loc[df['shield'] > 35] = 0
 |      >>> df
 |                  max_speed  shield
 |      cobra              30      10
 |      viper               0       0
 |      sidewinder          0       0
 |      
 |      **Getting values on a DataFrame with an index that has integer labels**
 |      
 |      Another example using integers for the index
 |      
 |      >>> df = pd.DataFrame([[1, 2], [4, 5], [7, 8]],
 |      ...      index=[7, 8, 9], columns=['max_speed', 'shield'])
 |      >>> df
 |         max_speed  shield
 |      7          1       2
 |      8          4       5
 |      9          7       8
 |      
 |      Slice with integer labels for rows. As mentioned above, note that both
 |      the start and stop of the slice are included.
 |      
 |      >>> df.loc[7:9]
 |         max_speed  shield
 |      7          1       2
 |      8          4       5
 |      9          7       8
 |      
 |      **Getting values with a MultiIndex**
 |      
 |      A number of examples using a DataFrame with a MultiIndex
 |      
 |      >>> tuples = [
 |      ...    ('cobra', 'mark i'), ('cobra', 'mark ii'),
 |      ...    ('sidewinder', 'mark i'), ('sidewinder', 'mark ii'),
 |      ...    ('viper', 'mark ii'), ('viper', 'mark iii')
 |      ... ]
 |      >>> index = pd.MultiIndex.from_tuples(tuples)
 |      >>> values = [[12, 2], [0, 4], [10, 20],
 |      ...         [1, 4], [7, 1], [16, 36]]
 |      >>> df = pd.DataFrame(values, columns=['max_speed', 'shield'], index=index)
 |      >>> df
 |                           max_speed  shield
 |      cobra      mark i           12       2
 |                 mark ii           0       4
 |      sidewinder mark i           10      20
 |                 mark ii           1       4
 |      viper      mark ii           7       1
 |                 mark iii         16      36
 |      
 |      Single label. Note this returns a DataFrame with a single index.
 |      
 |      >>> df.loc['cobra']
 |               max_speed  shield
 |      mark i          12       2
 |      mark ii          0       4
 |      
 |      Single index tuple. Note this returns a Series.
 |      
 |      >>> df.loc[('cobra', 'mark ii')]
 |      max_speed    0
 |      shield       4
 |      Name: (cobra, mark ii), dtype: int64
 |      
 |      Single label for row and column. Similar to passing in a tuple, this
 |      returns a Series.
 |      
 |      >>> df.loc['cobra', 'mark i']
 |      max_speed    12
 |      shield        2
 |      Name: (cobra, mark i), dtype: int64
 |      
 |      Single tuple. Note using ``[[]]`` returns a DataFrame.
 |      
 |      >>> df.loc[[('cobra', 'mark ii')]]
 |                     max_speed  shield
 |      cobra mark ii          0       4
 |      
 |      Single tuple for the index with a single label for the column
 |      
 |      >>> df.loc[('cobra', 'mark i'), 'shield']
 |      2
 |      
 |      Slice from index tuple to single label
 |      
 |      >>> df.loc[('cobra', 'mark i'):'viper']
 |                           max_speed  shield
 |      cobra      mark i           12       2
 |                 mark ii           0       4
 |      sidewinder mark i           10      20
 |                 mark ii           1       4
 |      viper      mark ii           7       1
 |                 mark iii         16      36
 |      
 |      Slice from index tuple to index tuple
 |      
 |      >>> df.loc[('cobra', 'mark i'):('viper', 'mark ii')]
 |                          max_speed  shield
 |      cobra      mark i          12       2
 |                 mark ii          0       4
 |      sidewinder mark i          10      20
 |                 mark ii          1       4
 |      viper      mark ii          7       1
 |      
 |      Raises
 |      ------
 |      KeyError:
 |          when any items are not found
 |  
 |  ndim
 |      Return an int representing the number of axes / array dimensions.
 |      
 |      Return 1 if Series. Otherwise return 2 if DataFrame.
 |      
 |      See Also
 |      --------
 |      ndarray.ndim
 |      
 |      Examples
 |      --------
 |      >>> s = pd.Series({'a': 1, 'b': 2, 'c': 3})
 |      >>> s.ndim
 |      1
 |      
 |      >>> df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
 |      >>> df.ndim
 |      2
 |  
 |  size
 |      Return an int representing the number of elements in this object.
 |      
 |      Return the number of rows if Series. Otherwise return the number of
 |      rows times number of columns if DataFrame.
 |      
 |      See Also
 |      --------
 |      ndarray.size
 |      
 |      Examples
 |      --------
 |      >>> s = pd.Series({'a': 1, 'b': 2, 'c': 3})
 |      >>> s.size
 |      3
 |      
 |      >>> df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
 |      >>> df.size
 |      4
 |  
 |  values
 |      Return a Numpy representation of the DataFrame.
 |      
 |      Only the values in the DataFrame will be returned, the axes labels
 |      will be removed.
 |      
 |      Returns
 |      -------
 |      numpy.ndarray
 |          The values of the DataFrame.
 |      
 |      Examples
 |      --------
 |      A DataFrame where all columns are the same type (e.g., int64) results
 |      in an array of the same type.
 |      
 |      >>> df = pd.DataFrame({'age':    [ 3,  29],
 |      ...                    'height': [94, 170],
 |      ...                    'weight': [31, 115]})
 |      >>> df
 |         age  height  weight
 |      0    3      94      31
 |      1   29     170     115
 |      >>> df.dtypes
 |      age       int64
 |      height    int64
 |      weight    int64
 |      dtype: object
 |      >>> df.values
 |      array([[  3,  94,  31],
 |             [ 29, 170, 115]], dtype=int64)
 |      
 |      A DataFrame with mixed type columns(e.g., str/object, int64, float32)
 |      results in an ndarray of the broadest type that accommodates these
 |      mixed types (e.g., object).
 |      
 |      >>> df2 = pd.DataFrame([('parrot',   24.0, 'second'),
 |      ...                     ('lion',     80.5, 1),
 |      ...                     ('monkey', np.nan, None)],
 |      ...                   columns=('name', 'max_speed', 'rank'))
 |      >>> df2.dtypes
 |      name          object
 |      max_speed    float64
 |      rank          object
 |      dtype: object
 |      >>> df2.values
 |      array([['parrot', 24.0, 'second'],
 |             ['lion', 80.5, 1],
 |             ['monkey', nan, None]], dtype=object)
 |      
 |      Notes
 |      -----
 |      The dtype will be a lower-common-denominator dtype (implicit
 |      upcasting); that is to say if the dtypes (even of numeric types)
 |      are mixed, the one that accommodates all will be chosen. Use this
 |      with care if you are not dealing with the blocks.
 |      
 |      e.g. If the dtypes are float16 and float32, dtype will be upcast to
 |      float32.  If dtypes are int32 and uint8, dtype will be upcast to
 |      int32. By :func:`numpy.find_common_type` convention, mixing int64
 |      and uint64 will result in a float64 dtype.
 |      
 |      See Also
 |      --------
 |      pandas.DataFrame.index : Retrievie the index labels
 |      pandas.DataFrame.columns : Retrieving the column names
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from pandas.core.base.PandasObject:
 |  
 |  __sizeof__(self)
 |      Generates the total memory usage for an object that returns
 |      either a value or Series of values
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from pandas.core.base.StringMixin:
 |  
 |  __bytes__(self)
 |      Return a string representation for a particular object.
 |      
 |      Invoked by bytes(obj) in py3 only.
 |      Yields a bytestring in both py2/py3.
 |  
 |  __repr__(self)
 |      Return a string representation for a particular object.
 |      
 |      Yields Bytestring in Py2, Unicode String in py3.
 |  
 |  __str__(self)
 |      Return a string representation for a particular Object
 |      
 |      Invoked by str(df) in both py2/py3.
 |      Yields Bytestring in Py2, Unicode String in py3.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors inherited from pandas.core.base.StringMixin:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
 |  
 |  ----------------------------------------------------------------------
 |  Methods inherited from pandas.core.accessor.DirNamesMixin:
 |  
 |  __dir__(self)
 |      Provide method name lookup and completion
 |      Only provide 'public' methods

```

## 2、软阈值计算

```pythonstart = datetime.datetime.now()
G_all=nx.Graph()
data_len=len(data)
data_list=data.index
#相关性系数
result=data.T.corr()
#运行时间
end = datetime.datetime.now()
print("时间：",end-start)```
*Output:*
```时间： 0:00:09.289201
```

```pythonnp.fill_diagonal(result.values, 0)
arr=abs(result)
arr_len=len(arr)
temp=np.zeros(arr_len)
for i in range(1,3):    
    temp=temp+(arr**i)/i
temp.head()```
*Output:*
```substanceBXH  MMT00000044  MMT00000046  MMT00000051  MMT00000080  MMT00000102  \
substanceBXH                                                                    
MMT00000044      0.000000     0.151256     0.127723     0.038499     0.174776   
MMT00000046      0.151256     0.000000     0.728334     0.049977     0.763454   
MMT00000051      0.127723     0.728334     0.000000     0.061974     0.096041   
MMT00000080      0.038499     0.049977     0.061974     0.000000     0.101928   
MMT00000102      0.174776     0.763454     0.096041     0.101928     0.000000   

substanceBXH  MMT00000149  MMT00000159  MMT00000207  MMT00000212  MMT00000241  \
substanceBXH                                                                    
MMT00000044      0.150469     0.039732     0.055749     0.040410     0.069559   
MMT00000046      0.296920     0.110256     0.194043     0.461836     0.083253   
MMT00000051      0.279626     0.719102     0.432563     0.549306     0.185569   
MMT00000080      0.118981     0.097212     0.597337     0.046650     0.838776   
MMT00000102      0.304612     0.449397     0.033212     0.099451     0.100090   

substanceBXH     ...       MMT00082753  MMT00082759  MMT00082798  MMT00082822  \
substanceBXH     ...                                                            
MMT00000044      ...          0.283788     0.259585     0.208980     0.162881   
MMT00000046      ...          0.386452     0.332012     0.578200     0.844994   
MMT00000051      ...          0.276479     0.040092     0.534604     0.705218   
MMT00000080      ...          0.137663     0.207288     0.012935     0.330876   
MMT00000102      ...          0.286503     0.301812     0.547826     0.265809   

substanceBXH  MMT00082829  MMT00082832  MMT00082850  MMT00082869  MMT00082877  \
substanceBXH                                                                    
MMT00000044      0.041975     0.034748     0.139239     0.148680     0.301681   
MMT00000046      0.221145     0.140971     0.628896     0.010682     0.325814   
MMT00000051      0.229944     0.510132     0.805841     0.070707     0.222062   
MMT00000080      1.164811     0.303138     0.067149     0.004718     0.152238   
MMT00000102      0.252107     0.174013     0.323302     0.288305     0.081244   

substanceBXH  MMT00082906  
substanceBXH               
MMT00000044      0.122667  
MMT00000046      0.075113  
MMT00000051      0.001092  
MMT00000080      0.699017  
MMT00000102      0.012353  

[5 rows x 2785 columns]```

```pythonre1=pd.DataFrame(columns=['beta','r2','meank'])
for j in range(1,12):
    result_i=np.float_power(temp,j)
    tt_0=np.sum(abs(result_i),axis=0)-1
    n=plt.hist(x = tt_0), # 指定绘图数据
    x=n[0][0]
    y=[]
    for i in range(len(n[0][1])-1):
        y.append((n[0][1][i]+n[0][1][i+1])/2)
    x=np.log10(x)
    y=np.log10(y)
    res=stats.linregress(x, y)
    r2=np.float_power(res.rvalue,2)
    k=tt_0.mean()
    re1=re1.append({'beta':j,'r2':r2,'meank':k},ignore_index=True)
re1```
*Output:*
```    beta        r2       meank
0    1.0  0.001579  683.970143
1    2.0  0.193164  297.709842
2    3.0  0.410875  173.086377
3    4.0  0.619393  120.827747
4    5.0  0.741184   96.126154
5    6.0  0.892932   84.551754
6    7.0  0.940257   80.597204
7    8.0  0.971655   82.057516
8    9.0  0.971214   88.225735
9   10.0  0.959321   99.251496
10  11.0  0.944997  115.930042```
```<Figure size 432x288 with 1 Axes>```

```pythonfor i in re1['r2']:
        if i>0.9:
            soft=re1[re1['r2']==i]['beta'].iloc[0]
            break
soft```
*Output:*
```7.0```

```pythonsns.color_palette("seismic")```
*Output:*
```[(0.0, 0.0, 0.6952941176470588),
 (0.1450980392156863, 0.1450980392156863, 1.0),
 (0.7098039215686275, 0.7098039215686275, 1.0),
 (1.0, 0.7098039215686274, 0.7098039215686274),
 (1.0, 0.14509803921568631, 0.14509803921568631),
 (0.7823529411764706, 0.0, 0.0)]```

```pythonsns.color_palette("mako")```
*Output:*
```[(0.18195582, 0.11955283, 0.23136943),
 (0.25307401, 0.23772973, 0.48316271),
 (0.21607792, 0.39736958, 0.61948028),
 (0.20344718, 0.56074869, 0.65649508),
 (0.25187832, 0.71827158, 0.67872193),
 (0.54578602, 0.8544913, 0.69848331)]```

```pythoncmap1[0]```
*Output:*
```(0.0, 0.0, 0.6952941176470588)```

```pythonhelp(sns.regplot)```
*Output:*
```Help on function regplot in module seaborn.regression:

regplot(x, y, data=None, x_estimator=None, x_bins=None, x_ci='ci', scatter=True, fit_reg=True, ci=95, n_boot=1000, units=None, order=1, logistic=False, lowess=False, robust=False, logx=False, x_partial=None, y_partial=None, truncate=False, dropna=True, x_jitter=None, y_jitter=None, label=None, color=None, marker='o', scatter_kws=None, line_kws=None, ax=None)
    Plot data and a linear regression model fit.
    
    There are a number of mutually exclusive options for estimating the
    regression model. See the :ref:`tutorial <regression_tutorial>` for more
    information.    
    
    Parameters
    ----------
    x, y: string, series, or vector array
        Input variables. If strings, these should correspond with column names
        in ``data``. When pandas objects are used, axes will be labeled with
        the series name.
    data : DataFrame
        Tidy ("long-form") dataframe where each column is a variable and each
        row is an observation.    
    x_estimator : callable that maps vector -> scalar, optional
        Apply this function to each unique value of ``x`` and plot the
        resulting estimate. This is useful when ``x`` is a discrete variable.
        If ``x_ci`` is given, this estimate will be bootstrapped and a
        confidence interval will be drawn.    
    x_bins : int or vector, optional
        Bin the ``x`` variable into discrete bins and then estimate the central
        tendency and a confidence interval. This binning only influences how
        the scatterplot is drawn; the regression is still fit to the original
        data.  This parameter is interpreted either as the number of
        evenly-sized (not necessary spaced) bins or the positions of the bin
        centers. When this parameter is used, it implies that the default of
        ``x_estimator`` is ``numpy.mean``.    
    x_ci : "ci", "sd", int in [0, 100] or None, optional
        Size of the confidence interval used when plotting a central tendency
        for discrete values of ``x``. If ``"ci"``, defer to the value of the
        ``ci`` parameter. If ``"sd"``, skip bootstrapping and show the
        standard deviation of the observations in each bin.    
    scatter : bool, optional
        If ``True``, draw a scatterplot with the underlying observations (or
        the ``x_estimator`` values).    
    fit_reg : bool, optional
        If ``True``, estimate and plot a regression model relating the ``x``
        and ``y`` variables.    
    ci : int in [0, 100] or None, optional
        Size of the confidence interval for the regression estimate. This will
        be drawn using translucent bands around the regression line. The
        confidence interval is estimated using a bootstrap; for large
        datasets, it may be advisable to avoid that computation by setting
        this parameter to None.    
    n_boot : int, optional
        Number of bootstrap resamples used to estimate the ``ci``. The default
        value attempts to balance time and stability; you may want to increase
        this value for "final" versions of plots.    
    units : variable name in ``data``, optional
        If the ``x`` and ``y`` observations are nested within sampling units,
        those can be specified here. This will be taken into account when
        computing the confidence intervals by performing a multilevel bootstrap
        that resamples both units and observations (within unit). This does not
        otherwise influence how the regression is estimated or drawn.    
    order : int, optional
        If ``order`` is greater than 1, use ``numpy.polyfit`` to estimate a
        polynomial regression.    
    logistic : bool, optional
        If ``True``, assume that ``y`` is a binary variable and use
        ``statsmodels`` to estimate a logistic regression model. Note that this
        is substantially more computationally intensive than linear regression,
        so you may wish to decrease the number of bootstrap resamples
        (``n_boot``) or set ``ci`` to None.    
    lowess : bool, optional
        If ``True``, use ``statsmodels`` to estimate a nonparametric lowess
        model (locally weighted linear regression). Note that confidence
        intervals cannot currently be drawn for this kind of model.    
    robust : bool, optional
        If ``True``, use ``statsmodels`` to estimate a robust regression. This
        will de-weight outliers. Note that this is substantially more
        computationally intensive than standard linear regression, so you may
        wish to decrease the number of bootstrap resamples (``n_boot``) or set
        ``ci`` to None.    
    logx : bool, optional
        If ``True``, estimate a linear regression of the form y ~ log(x), but
        plot the scatterplot and regression model in the input space. Note that
        ``x`` must be positive for this to work.    
    {x,y}_partial : strings in ``data`` or matrices
        Confounding variables to regress out of the ``x`` or ``y`` variables
        before plotting.    
    truncate : bool, optional
        By default, the regression line is drawn to fill the x axis limits
        after the scatterplot is drawn. If ``truncate`` is ``True``, it will
        instead by bounded by the data limits.    
    {x,y}_jitter : floats, optional
        Add uniform random noise of this size to either the ``x`` or ``y``
        variables. The noise is added to a copy of the data after fitting the
        regression, and only influences the look of the scatterplot. This can
        be helpful when plotting variables that take discrete values.    
    label : string
        Label to apply to ether the scatterplot or regression line (if
        ``scatter`` is ``False``) for use in a legend.
    color : matplotlib color
        Color to apply to all plot elements; will be superseded by colors
        passed in ``scatter_kws`` or ``line_kws``.
    marker : matplotlib marker code
        Marker to use for the scatterplot glyphs.
    {scatter,line}_kws : dictionaries
        Additional keyword arguments to pass to ``plt.scatter`` and
        ``plt.plot``.    
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    
    Returns
    -------
    ax : matplotlib Axes
        The Axes object containing the plot.
    
    See Also
    --------
    lmplot : Combine :func:`regplot` and :class:`FacetGrid` to plot multiple
             linear relationships in a dataset.
    jointplot : Combine :func:`regplot` and :class:`JointGrid` (when used with
                ``kind="reg"``).
    pairplot : Combine :func:`regplot` and :class:`PairGrid` (when used with
               ``kind="reg"``).
    residplot : Plot the residuals of a linear regression model.
    
    Notes
    -----
    
    The :func:`regplot` and :func:`lmplot` functions are closely related, but
    the former is an axes-level function while the latter is a figure-level
    function that combines :func:`regplot` and :class:`FacetGrid`.    
    
    
    It's also easy to combine combine :func:`regplot` and :class:`JointGrid` or
    :class:`PairGrid` through the :func:`jointplot` and :func:`pairplot`
    functions, although these do not directly accept all of :func:`regplot`'s
    parameters.
    
    Examples
    --------
    
    Plot the relationship between two variables in a DataFrame:
    
    .. plot::
        :context: close-figs
    
        >>> import seaborn as sns; sns.set(color_codes=True)
        >>> tips = sns.load_dataset("tips")
        >>> ax = sns.regplot(x="total_bill", y="tip", data=tips)
    
    Plot with two variables defined as numpy arrays; use a different color:
    
    .. plot::
        :context: close-figs
    
        >>> import numpy as np; np.random.seed(8)
        >>> mean, cov = [4, 6], [(1.5, .7), (.7, 1)]
        >>> x, y = np.random.multivariate_normal(mean, cov, 80).T
        >>> ax = sns.regplot(x=x, y=y, color="g")
    
    Plot with two variables defined as pandas Series; use a different marker:
    
    .. plot::
        :context: close-figs
    
        >>> import pandas as pd
        >>> x, y = pd.Series(x, name="x_var"), pd.Series(y, name="y_var")
        >>> ax = sns.regplot(x=x, y=y, marker="+")
    
    Use a 68% confidence interval, which corresponds with the standard error
    of the estimate:
    
    .. plot::
        :context: close-figs
    
        >>> ax = sns.regplot(x=x, y=y, ci=68)
    
    Plot with a discrete ``x`` variable and add some jitter:
    
    .. plot::
        :context: close-figs
    
        >>> ax = sns.regplot(x="size", y="total_bill", data=tips, x_jitter=.1)
    
    Plot with a discrete ``x`` variable showing means and confidence intervals
    for unique values:
    
    .. plot::
        :context: close-figs
    
        >>> ax = sns.regplot(x="size", y="total_bill", data=tips,
        ...                  x_estimator=np.mean)
    
    Plot with a continuous variable divided into discrete bins:
    
    .. plot::
        :context: close-figs
    
        >>> ax = sns.regplot(x=x, y=y, x_bins=4)
    
    Fit a higher-order polynomial regression and truncate the model prediction:
    
    .. plot::
        :context: close-figs
    
        >>> ans = sns.load_dataset("anscombe")
        >>> ax = sns.regplot(x="x", y="y", data=ans.loc[ans.dataset == "II"],
        ...                  scatter_kws={"s": 80},
        ...                  order=2, ci=None, truncate=True)
    
    Fit a robust regression and don't plot a confidence interval:
    
    .. plot::
        :context: close-figs
    
        >>> ax = sns.regplot(x="x", y="y", data=ans.loc[ans.dataset == "III"],
        ...                  scatter_kws={"s": 80},
        ...                  robust=True, ci=None)
    
    Fit a logistic regression; jitter the y variable and use fewer bootstrap
    iterations:
    
    .. plot::
        :context: close-figs
    
        >>> tips["big_tip"] = (tips.tip / tips.total_bill) > .175
        >>> ax = sns.regplot(x="total_bill", y="big_tip", data=tips,
        ...                  logistic=True, n_boot=500, y_jitter=.03)
    
    Fit the regression model using log(x) and truncate the model prediction:
    
    .. plot::
        :context: close-figs
    
        >>> ax = sns.regplot(x="size", y="total_bill", data=tips,
        ...                  x_estimator=np.mean, logx=True, truncate=True)

```

```python
my_dpi=300
fig=plt.figure(figsize=(2000/my_dpi, 1000/my_dpi), dpi=my_dpi)
cc='seismic'
cmap = sns.color_palette(cc)
grid = plt.GridSpec(1, 4, wspace=1, hspace=0.1)
#fig, (ax0, ax1) = plt.subplots(2, 1)
plt.subplot(grid[0,0:2])
p1=sns.regplot(x=re1["beta"], y=re1['r2'], fit_reg=False, marker="o", color="#42B520")
p1.axhline(y=0.9,ls=":",c=cmap[5])

#fig, (ax0, ax1) = plt.subplots(2, 1)
plt.subplot(grid[0,2:4])
p1=sns.regplot(x=re1["beta"], y=re1['meank'], fit_reg=False, marker="o", color=cmap[4])

plt.show()```
*Output:*
```<Figure size 2000x1000 with 2 Axes>```

```pythontt_0=np.sum(abs(np.float_power(temp,6)),axis=0)-1
plt.hist(x = tt_0), # 指定绘图数据```
*Output:*
```((array([1211.,  695.,  397.,  215.,  105.,   60.,   34.,   37.,   17.,
           14.]),
  array([ -0.4909389 ,  47.29907759,  95.08909407, 142.87911055,
         190.66912704, 238.45914352, 286.24916001, 334.03917649,
         381.82919297, 429.61920946, 477.40922594]),
  <a list of 10 Patch objects>),)```
```<Figure size 432x288 with 1 Axes>```

```pythontest```
*Output:*
```substanceBXH   MMT00000044   MMT00000046   MMT00000051   MMT00000080  \
substanceBXH                                                           
MMT00000044   1.000000e+00  1.197504e-05  4.341150e-06  3.256149e-09   
MMT00000046   1.197504e-05  1.000000e+00  1.492743e-01  1.558217e-08   
MMT00000051   4.341150e-06  1.492743e-01  1.000000e+00  5.665705e-08   
MMT00000080   3.256149e-09  1.558217e-08  5.665705e-08  1.000000e+00   
MMT00000102   2.850266e-05  1.980153e-01  7.847442e-07  1.121383e-06   
MMT00000149   1.160594e-05  6.852361e-04  4.780371e-04  2.837109e-06   
MMT00000159   3.934096e-09  1.796414e-06  1.382745e-01  8.439377e-07   
MMT00000207   3.001930e-08  5.338055e-05  6.550819e-03  4.542708e-02   
MMT00000212   4.354126e-09  9.703476e-03  2.747182e-02  1.030636e-08   
MMT00000241   1.132742e-07  3.329632e-07  4.083473e-05  3.482379e-01   
MMT00000268   3.206218e-06  1.131007e-08  1.502435e-09  6.318836e-01   
MMT00000283   9.491159e-07  2.011449e-05  1.044093e-09  2.269627e-04   
MMT00000334   2.630881e-07  3.440183e-01  1.379453e-01  2.749600e-04   
MMT00000373   4.130444e-05  1.315786e-07  4.933695e-07  2.974199e-05   
MMT00000401   2.347990e-11  2.657256e-11  6.519745e-04  1.754935e-05   
MMT00000418   2.234372e-20  3.514663e-07  3.062768e-05  1.251110e-04   
MMT00000464   4.862072e-07  3.552448e-01  1.310801e-01  4.770903e-08   
MMT00000517   8.234659e-06  1.407421e-02  4.619389e-02  3.321338e-04   
MMT00000549   4.169165e-08  2.210103e-07  2.634914e-04  1.302562e-03   
MMT00000602   5.598813e-06  1.707741e-01  6.871978e-03  1.618902e-04   
MMT00000608   2.327595e-05  6.272342e-04  1.180358e-06  4.811260e-07   
MMT00000701   6.680642e-06  1.032428e-03  3.833191e-04  2.138193e-06   
MMT00000713   8.481782e-07  1.296452e-07  1.360947e-13  2.086756e-10   
MMT00000743   4.501496e-10  1.344219e-07  4.036721e-04  6.532767e-07   
MMT00000793   3.410415e-05  2.105562e-04  2.094335e-01  2.093537e-03   
MMT00000840   5.661111e-14  5.647709e-04  5.924822e-06  6.818126e-05   
MMT00000887   8.118060e-06  1.445379e-03  9.778844e-06  1.147039e-06   
MMT00000963   2.716186e-04  1.662243e-02  2.063618e-04  3.646373e-10   
MMT00000988   9.791968e-08  5.847439e-06  5.658489e-06  1.208775e-09   
MMT00000996   2.724761e-08  3.809444e-02  1.587004e-02  4.171269e-05   
...                    ...           ...           ...           ...   
MMT00082101   1.632469e-09  2.569959e-01  5.885682e-02  1.201304e-05   
MMT00082110   5.168364e-04  1.214540e-09  2.135647e-07  4.162376e-06   
MMT00082126   6.501050e-06  1.781455e-04  1.898983e-03  1.237950e-04   
MMT00082164   1.879628e-05  6.080309e-01  3.513716e-01  2.118398e-03   
MMT00082181   9.150493e-06  1.008763e-04  2.017434e-02  2.105940e-04   
MMT00082259   6.280820e-05  1.253572e+00  1.468325e-01  7.674771e-04   
MMT00082420   3.968635e-11  5.470475e-07  5.073998e-04  6.816866e-06   
MMT00082428   3.221176e-07  2.335142e-07  1.191083e-04  3.338304e-08   
MMT00082445   7.730947e-07  4.936872e-06  1.061072e-09  7.448434e-05   
MMT00082461   1.886821e-07  1.626089e-01  2.907498e-01  1.124044e-06   
MMT00082551   6.501764e-09  2.634476e-04  7.365689e-02  2.233470e-06   
MMT00082577   2.230942e-08  2.432304e-03  2.393712e-02  2.050249e-03   
MMT00082579   2.689731e-05  1.296375e+00  4.864903e-01  1.427994e-08   
MMT00082585   9.283044e-04  1.382598e-01  1.197970e-02  1.330376e-07   
MMT00082592   2.598686e-05  3.084140e-01  1.540391e-01  3.533306e-06   
MMT00082622   2.024555e-07  2.514624e-03  8.431933e-06  2.147554e-02   
MMT00082651   4.493824e-04  1.119592e-04  6.814406e-04  1.178374e-03   
MMT00082663   6.641392e-05  6.691464e-01  2.589723e-01  6.130081e-06   
MMT00082677   9.286917e-10  4.539510e-04  2.817543e-02  3.131211e-03   
MMT00082712   2.304455e-06  4.437369e-06  3.473246e-02  4.185850e-05   
MMT00082753   5.223542e-04  3.330974e-03  4.466520e-04  6.806124e-06   
MMT00082759   3.059698e-04  1.339449e-03  4.152634e-09  7.933190e-05   
MMT00082798   8.329766e-05  3.736526e-02  2.334498e-02  4.683582e-12   
MMT00082822   1.867339e-05  3.640183e-01  1.230099e-01  1.312171e-03   
MMT00082829   5.469334e-09  1.169656e-04  1.478197e-04  2.497652e+00   
MMT00082832   1.760377e-09  7.848353e-06  1.762360e-02  7.759591e-04   
MMT00082850   7.287326e-06  6.186918e-02  2.738400e-01  9.167311e-08   
MMT00082869   1.080208e-05  1.486003e-12  1.249611e-07  1.102555e-14   
MMT00082877   7.538520e-04  1.196235e-03  1.199083e-04  1.244891e-05   
MMT00082906   3.406891e-06  1.795923e-07  1.692309e-18  1.166613e-01   

substanceBXH   MMT00000102   MMT00000149   MMT00000159   MMT00000207  \
substanceBXH                                                           
MMT00000044   2.850266e-05  1.160594e-05  3.934096e-09  3.001930e-08   
MMT00000046   1.980153e-01  6.852361e-04  1.796414e-06  5.338055e-05   
MMT00000051   7.847442e-07  4.780371e-04  1.382745e-01  6.550819e-03   
MMT00000080   1.121383e-06  2.837109e-06  8.439377e-07  4.542708e-02   
MMT00000102   1.000000e+00  7.988876e-04  8.237175e-03  1.342171e-09   
MMT00000149   7.988876e-04  1.000000e+00  3.659330e-06  6.025893e-02   
MMT00000159   8.237175e-03  3.659330e-06  1.000000e+00  1.308707e-02   
MMT00000207   1.342171e-09  6.025893e-02  1.308707e-02  1.000000e+00   
MMT00000212   9.674992e-07  8.914365e-05  4.457971e-03  1.091092e-03   
MMT00000241   1.005440e-06  8.861248e-03  7.788086e-05  9.313010e-02   
MMT00000268   1.850107e-08  6.964226e-05  5.544693e-08  4.629398e-03   
MMT00000283   1.284202e-06  1.634957e-01  4.508556e-06  1.969352e-03   
MMT00000334   2.528099e-03  7.792270e-03  2.149511e-04  8.385169e-03   
MMT00000373   2.813143e-04  4.296860e-06  1.456766e-04  2.740717e-07   
MMT00000401   4.281944e-05  4.453617e-06  1.080934e-02  6.555697e-05   
MMT00000418   1.411375e-07  1.434849e-01  7.759098e-06  2.754895e-03   
MMT00000464   1.678568e-03  4.536184e-04  1.689350e-03  1.918417e-03   
MMT00000517   1.539860e-04  9.873679e-03  8.629538e-04  4.116094e-02   
MMT00000549   1.806666e-06  3.416601e-06  2.566497e-03  3.216206e-03   
MMT00000602   7.324179e-04  1.515238e-06  2.031745e-06  1.350798e-06   
MMT00000608   1.670028e-03  8.874526e-10  1.285467e-05  9.408082e-12   
MMT00000701   3.082795e-05  1.783162e-03  4.370752e-13  2.355403e-07   
MMT00000713   5.554289e-06  3.588038e-11  7.640847e-06  1.583919e-05   
MMT00000743   1.102813e-04  6.044102e-04  2.821708e-03  1.616809e-03   
MMT00000793   1.992957e-10  6.551391e-03  2.083653e-01  9.861808e-02   
MMT00000840   1.079141e-02  8.286070e-10  9.856421e-03  4.307060e-03   
MMT00000887   1.825844e-01  1.608086e-04  4.480274e-03  6.110802e-06   
MMT00000963   5.256529e-03  5.184465e-05  6.053996e-07  5.129655e-11   
MMT00000988   5.894979e-04  8.738387e-07  8.852542e-03  4.901389e-06   
MMT00000996   5.884650e-07  4.012047e-03  1.887000e-03  3.716125e-02   
...                    ...           ...           ...           ...   
MMT00082101   3.050283e-02  5.767326e-03  4.773243e-05  5.071081e-04   
MMT00082110   1.788061e-08  1.568401e-09  3.503041e-05  7.285250e-08   
MMT00082126   2.173230e-06  2.473934e-01  6.048031e-04  3.263439e-02   
MMT00082164   1.855317e-03  2.917246e-04  4.253931e-03  4.349320e-06   
MMT00082181   5.624263e-07  7.977066e-08  1.174612e-03  9.699006e-06   
MMT00082259   3.515062e-03  6.170286e-03  3.530593e-03  4.993699e-05   
MMT00082420   2.603235e-03  5.107258e-11  1.209558e-02  2.099414e-03   
MMT00082428   4.628707e-03  9.033662e-07  1.244371e-02  1.160422e-05   
MMT00082445   3.177401e-03  8.917804e-08  3.116602e-03  4.577868e-09   
MMT00082461   7.533443e-07  9.042689e-03  4.904726e-02  7.669263e-02   
MMT00082551   4.712344e-05  7.620381e-15  8.860336e-02  1.223814e-03   
MMT00082577   6.091210e-19  1.132535e-01  4.024124e-03  7.463448e-01   
MMT00082579   2.415178e-02  4.982600e-02  1.304134e-04  2.348113e-03   
MMT00082585   1.214721e-02  6.179322e-05  1.876841e-07  6.678118e-06   
MMT00082592   4.932935e-03  4.116434e-02  1.439763e-04  3.708175e-06   
MMT00082622   1.549872e-02  2.806117e-01  1.640206e-05  8.528350e-03   
MMT00082651   1.149918e-07  4.671528e-02  3.692595e-07  2.164902e-05   
MMT00082663   2.342246e-05  3.738869e-02  4.904988e-03  1.109477e-02   
MMT00082677   2.978438e-08  1.973139e-05  6.346309e-03  1.832332e-02   
MMT00082712   1.943379e-04  3.103831e-04  6.400580e-02  2.269388e-04   
MMT00082753   5.530631e-04  1.084290e-03  8.335968e-08  2.419597e-08   
MMT00082759   7.558232e-04  5.028469e-20  8.053848e-07  1.580722e-05   
MMT00082798   2.703055e-02  2.407982e-02  1.165694e-08  3.966177e-05   
MMT00082822   3.527133e-04  2.022625e-14  5.364075e-03  7.835153e-08   
MMT00082829   2.567520e-04  6.039651e-03  5.549562e-05  2.292926e-01   
MMT00082832   2.776416e-05  2.552469e-07  3.490221e-02  2.691238e-03   
MMT00082850   1.141964e-03  8.921925e-06  8.247591e-05  3.535673e-05   
MMT00082869   5.742685e-04  6.106139e-02  1.183968e-05  2.437774e-08   
MMT00082877   2.875780e-07  5.063917e-02  2.008927e-08  1.029818e-02   
MMT00082906   3.554065e-12  8.391253e-04  2.011392e-07  1.362075e-01   

substanceBXH   MMT00000212   MMT00000241      ...        MMT00082753  \
substanceBXH                                  ...                      
MMT00000044   4.354126e-09  1.132742e-07      ...       5.223542e-04   
MMT00000046   9.703476e-03  3.329632e-07      ...       3.330974e-03   
MMT00000051   2.747182e-02  4.083473e-05      ...       4.466520e-04   
MMT00000080   1.030636e-08  3.482379e-01      ...       6.806124e-06   
MMT00000102   9.674992e-07  1.005440e-06      ...       5.530631e-04   
MMT00000149   8.914365e-05  8.861248e-03      ...       1.084290e-03   
MMT00000159   4.457971e-03  7.788086e-05      ...       8.335968e-08   
MMT00000207   1.091092e-03  9.313010e-02      ...       2.419597e-08   
MMT00000212   1.000000e+00  2.045567e-04      ...       7.513430e-05   
MMT00000241   2.045567e-04  1.000000e+00      ...       5.128607e-06   
MMT00000268   8.429458e-07  2.771651e+00      ...       1.545941e-06   
MMT00000283   2.334301e-06  1.602512e-14      ...       1.161035e-05   
MMT00000334   3.009126e-02  1.669498e-02      ...       1.207354e-02   
MMT00000373   6.898117e-04  8.793581e-07      ...       8.815757e-08   
MMT00000401   4.594919e-06  9.283754e-06      ...       1.191590e-07   
MMT00000418   4.156965e-07  3.691830e-03      ...       9.442362e-07   
MMT00000464   2.180824e-02  3.302567e-07      ...       1.164424e-03   
MMT00000517   2.164666e-02  7.510875e-04      ...       2.873312e-03   
MMT00000549   4.101594e-03  2.470266e-03      ...       1.692828e-12   
MMT00000602   2.492721e-03  2.899936e-04      ...       4.188470e-03   
MMT00000608   6.128014e-06  2.028801e-05      ...       7.559177e-04   
MMT00000701   1.147835e-04  4.195466e-07      ...       6.496845e-02   
MMT00000713   4.264350e-08  6.022708e-15      ...       8.861098e-08   
MMT00000743   5.466572e-04  2.799786e-03      ...       3.805865e-09   
MMT00000793   6.658256e-03  2.707770e-02      ...       8.117825e-04   
MMT00000840   5.084727e-05  1.125927e-04      ...       5.444215e-09   
MMT00000887   1.680516e-05  2.270564e-08      ...       1.042094e-06   
MMT00000963   3.918678e-06  4.620232e-06      ...       5.157728e-04   
MMT00000988   1.274134e-03  6.840251e-08      ...       7.679692e-06   
MMT00000996   8.225533e-03  2.610542e-04      ...       3.291872e-05   
...                    ...           ...      ...                ...   
MMT00082101   3.649090e-02  1.450694e-03      ...       4.544629e-03   
MMT00082110   5.474205e-10  1.935183e-08      ...       7.480623e-07   
MMT00082126   7.948415e-05  1.027922e-03      ...       5.669973e-04   
MMT00082164   6.833385e-02  1.016683e-05      ...       1.568561e-03   
MMT00082181   1.068352e-04  3.027398e-04      ...       2.975320e-08   
MMT00082259   1.954935e-02  4.761191e-09      ...       3.115902e-03   
MMT00082420   4.248240e-04  1.372631e-13      ...       7.290044e-10   
MMT00082428   2.638875e-04  1.421187e-05      ...       2.311399e-06   
MMT00082445   9.347011e-05  1.408855e-03      ...       1.682401e-05   
MMT00082461   1.202066e-01  8.097443e-04      ...       3.610423e-03   
MMT00082551   2.533429e-03  6.829315e-04      ...       1.580215e-05   
MMT00082577   1.393435e-03  4.997952e-03      ...       8.478964e-06   
MMT00082579   6.051916e-02  4.973794e-05      ...       2.495276e-02   
MMT00082585   2.263373e-03  2.471552e-07      ...       1.522245e-02   
MMT00082592   6.051460e-04  1.181690e-10      ...       1.183277e-03   
MMT00082622   7.096740e-12  5.703954e-02      ...       6.193839e-06   
MMT00082651   1.043666e-13  1.438629e-05      ...       1.586001e-03   
MMT00082663   1.117301e-02  6.355644e-14      ...       4.352152e-03   
MMT00082677   3.433715e-02  2.887781e-02      ...       4.018311e-06   
MMT00082712   6.275460e-03  1.973649e-06      ...       1.529175e-14   
MMT00082753   7.513430e-05  5.128607e-06      ...       1.000000e+00   
MMT00082759   8.469747e-10  6.273649e-05      ...       1.681465e-07   
MMT00082798   6.605621e-06  2.371423e-10      ...       4.352261e-04   
MMT00082822   5.180256e-02  3.122564e-05      ...       3.123091e-04   
MMT00082829   2.949492e-05  2.492592e-01      ...       1.122071e-12   
MMT00082832   6.947933e-04  9.187112e-07      ...       2.234119e-08   
MMT00082850   4.808289e-05  1.301599e-10      ...       5.115783e-05   
MMT00082869   6.027643e-04  1.080550e-06      ...       1.203741e-13   
MMT00082877   4.064906e-07  5.917055e-04      ...       5.574431e-05   
MMT00082906   2.771614e-07  3.077123e-03      ...       2.874204e-04   

substanceBXH   MMT00082759   MMT00082798   MMT00082822   MMT00082829  \
substanceBXH                                                           
MMT00000044   3.059698e-04  8.329766e-05  1.867339e-05  5.469334e-09   
MMT00000046   1.339449e-03  3.736526e-02  3.640183e-01  1.169656e-04   
MMT00000051   4.152634e-09  2.334498e-02  1.230099e-01  1.478197e-04   
MMT00000080   7.933190e-05  4.683582e-12  1.312171e-03  2.497652e+00   
MMT00000102   7.558232e-04  2.703055e-02  3.527133e-04  2.567520e-04   
MMT00000149   5.028469e-20  2.407982e-02  2.022625e-14  6.039651e-03   
MMT00000159   8.053848e-07  1.165694e-08  5.364075e-03  5.549562e-05   
MMT00000207   1.580722e-05  3.966177e-05  7.835153e-08  2.292926e-01   
MMT00000212   8.469747e-10  6.605621e-06  5.180256e-02  2.949492e-05   
MMT00000241   6.273649e-05  2.371423e-10  3.122564e-05  2.492592e-01   
MMT00000268   1.157509e-04  7.015701e-06  2.608728e-03  1.068502e-01   
MMT00000283   9.189850e-06  4.399116e-05  1.057982e-04  8.497479e-10   
MMT00000334   1.383009e-07  2.712031e-03  5.035711e-03  1.185087e-02   
MMT00000373   2.400990e-05  4.126672e-10  2.960808e-11  9.640134e-07   
MMT00000401   2.664783e-06  3.238031e-04  4.690562e-04  1.040767e-07   
MMT00000418   3.895280e-08  9.958259e-04  3.223302e-02  3.229739e-03   
MMT00000464   2.121656e-09  4.259248e-02  1.383615e-01  2.512925e-03   
MMT00000517   1.736726e-08  4.422692e-05  7.261856e-03  6.981217e-03   
MMT00000549   9.575914e-05  2.029772e-06  7.230211e-04  2.571495e-04   
MMT00000602   2.552409e-05  9.648867e-09  7.826084e-02  1.542744e-08   
MMT00000608   2.761253e-05  7.319839e-06  2.136679e-06  1.912661e-06   
MMT00000701   1.585160e-12  6.253606e-06  2.022291e-07  2.783350e-12   
MMT00000713   1.914103e-04  7.064698e-06  4.580934e-10  6.251029e-08   
MMT00000743   3.644338e-17  4.778122e-05  6.358967e-06  1.007573e-08   
MMT00000793   2.522444e-09  7.941749e-05  2.179849e-05  8.376866e-03   
MMT00000840   1.565477e-06  2.057243e-08  1.520239e-06  5.435769e-06   
MMT00000887   9.993056e-04  1.632688e-03  3.969470e-08  1.605212e-11   
MMT00000963   1.761801e-02  1.622480e-03  9.962700e-05  2.324009e-07   
MMT00000988   6.496244e-05  8.702054e-09  4.347216e-09  5.045560e-06   
MMT00000996   1.741160e-05  1.236681e-04  1.564908e-02  4.841508e-03   
...                    ...           ...           ...           ...   
MMT00082101   4.393964e-06  4.247792e-04  3.522087e-02  5.211724e-04   
MMT00082110   1.690211e-08  2.929405e-11  1.263524e-11  6.973103e-06   
MMT00082126   1.115825e-05  4.268799e-03  4.834141e-17  6.850922e-03   
MMT00082164   8.079399e-05  5.644218e-03  2.553266e+00  4.733836e-07   
MMT00082181   1.921933e-05  2.086614e-02  3.496830e-03  2.800827e-10   
MMT00082259   5.333831e-04  3.863312e-02  9.096065e-01  1.999457e-12   
MMT00082420   1.877343e-04  2.195278e-07  2.420470e-08  1.800191e-05   
MMT00082428   3.041892e-04  3.787141e-06  1.209182e-08  3.614272e-07   
MMT00082445   2.002141e-03  2.025112e-08  3.360474e-08  1.359508e-07   
MMT00082461   5.690015e-06  2.497230e-04  1.267624e-01  1.126729e-03   
MMT00082551   1.664996e-10  4.632959e-07  1.188794e-02  1.677930e-05   
MMT00082577   2.549666e-06  2.008042e-03  1.702018e-06  1.584526e-01   
MMT00082579   2.320018e-04  1.419185e-01  3.336523e-01  1.890448e-04   
MMT00082585   8.264894e-06  1.599204e-06  6.670448e-03  2.258722e-05   
MMT00082592   2.228028e-04  2.243935e-01  7.549574e-03  2.101132e-06   
MMT00082622   1.443246e-06  7.894867e-03  4.322414e-04  1.803313e-01   
MMT00082651   1.227082e-04  1.532498e-01  2.891435e-09  9.219567e-14   
MMT00082663   1.797688e-05  1.455159e-01  8.094472e-02  1.876318e-04   
MMT00082677   1.294779e-07  4.298702e-06  7.784637e-04  1.364526e-03   
MMT00082712   7.305365e-04  2.646345e-06  2.969784e-02  9.260870e-07   
MMT00082753   1.681465e-07  4.352261e-04  3.123091e-04  1.122071e-12   
MMT00082759   1.000000e+00  1.372639e-03  2.169374e-04  1.234088e-05   
MMT00082798   1.372639e-03  1.000000e+00  5.796678e-04  1.367737e-03   
MMT00082822   2.169374e-04  5.796678e-04  1.000000e+00  4.705280e-06   
MMT00082829   1.234088e-05  1.367737e-03  4.705280e-06  1.000000e+00   
MMT00082832   1.370091e-03  1.315196e-06  1.123074e-05  3.227641e-03   
MMT00082850   2.215884e-04  1.149304e+00  1.407546e-02  6.135732e-04   
MMT00082869   7.991345e-06  2.788908e-02  9.438444e-05  4.909215e-06   
MMT00082877   6.695043e-09  5.653087e-03  2.004508e-08  2.129642e-03   
MMT00082906   2.221897e-06  1.069492e-06  1.583379e-05  1.462919e-01   

substanceBXH   MMT00082832   MMT00082850   MMT00082869   MMT00082877  \
substanceBXH                                                           
MMT00000044   1.760377e-09  7.287326e-06  1.080208e-05  7.538520e-04   
MMT00000046   7.848353e-06  6.186918e-02  1.486003e-12  1.196235e-03   
MMT00000051   1.762360e-02  2.738400e-01  1.249611e-07  1.199083e-04   
MMT00000080   7.759591e-04  9.167311e-08  1.102555e-14  1.244891e-05   
MMT00000102   2.776416e-05  1.141964e-03  5.742685e-04  2.875780e-07   
MMT00000149   2.552469e-07  8.921925e-06  6.106139e-02  5.063917e-02   
MMT00000159   3.490221e-02  8.247591e-05  1.183968e-05  2.008927e-08   
MMT00000207   2.691238e-03  3.535673e-05  2.437774e-08  1.029818e-02   
MMT00000212   6.947933e-04  4.808289e-05  6.027643e-04  4.064906e-07   
MMT00000241   9.187112e-07  1.301599e-10  1.080550e-06  5.917055e-04   
MMT00000268   4.924194e-11  1.876122e-06  2.086066e-09  5.568058e-05   
MMT00000283   2.701470e-06  8.252280e-05  9.813822e-03  1.661415e-03   
MMT00000334   1.815988e-03  1.917814e-02  8.813731e-06  2.073215e-03   
MMT00000373   2.885540e-04  6.942717e-10  2.924892e-09  4.314991e-10   
MMT00000401   6.748131e-03  2.428712e-04  2.589087e-10  2.376860e-12   
MMT00000418   4.450421e-07  2.544948e-08  5.415050e-02  2.615130e-02   
MMT00000464   6.418535e-02  6.907624e-02  6.049190e-07  1.543352e-04   
MMT00000517   3.434292e-04  6.026788e-05  1.314384e-05  8.126470e-04   
MMT00000549   1.968368e-03  4.258944e-08  1.477194e-02  6.858324e-09   
MMT00000602   6.128027e-07  2.954122e-07  6.567852e-04  2.309068e-06   
MMT00000608   1.452240e-06  3.147574e-06  2.669293e-07  4.367271e-06   
MMT00000701   2.202296e-09  7.308502e-09  2.628887e-07  1.040085e-04   
MMT00000713   9.821334e-08  7.266192e-06  2.928110e-09  5.678956e-09   
MMT00000743   1.479739e-05  8.804604e-06  1.693622e-07  2.427351e-04   
MMT00000793   1.441092e-02  2.757416e-03  1.160544e-07  2.825614e-03   
MMT00000840   6.245231e-05  2.639399e-13  1.415391e-03  1.296108e-09   
MMT00000887   4.044265e-02  2.744024e-09  1.648141e-01  8.841053e-06   
MMT00000963   7.469335e-04  1.142065e-04  1.629769e-03  1.066767e-03   
MMT00000988   9.911793e-03  6.783529e-08  3.176910e-04  1.400638e-07   
MMT00000996   3.443589e-02  6.170835e-04  2.045449e-05  2.661750e-03   
...                    ...           ...           ...           ...   
MMT00082101   4.402602e-04  3.561777e-03  1.460012e-05  7.454124e-07   
MMT00082110   4.753989e-08  9.056133e-08  5.791211e-13  2.771459e-06   
MMT00082126   1.149368e-05  1.354001e-05  1.061532e-04  1.538684e-02   
MMT00082164   3.729930e-05  2.990597e-02  2.513575e-09  1.946595e-07   
MMT00082181   9.439104e-04  2.031848e-02  1.116600e-09  1.559387e-04   
MMT00082259   3.135926e-04  3.517273e-02  3.652989e-06  8.278834e-04   
MMT00082420   1.141817e+00  2.727809e-07  3.013196e-04  1.427762e-06   
MMT00082428   1.105388e+00  1.700722e-10  1.138696e-03  8.847564e-11   
MMT00082445   2.362702e-03  3.651710e-08  5.913280e-06  3.260577e-07   
MMT00082461   1.216657e-01  7.560028e-03  3.923497e-04  1.537145e-03   
MMT00082551   2.226462e-02  3.253282e-03  5.759029e-04  4.881561e-07   
MMT00082577   1.357006e-02  1.010947e-03  3.576590e-09  1.599913e-01   
MMT00082579   1.244693e-04  1.536536e-01  3.802010e-08  9.158608e-04   
MMT00082585   6.234754e-07  1.612617e-05  9.812385e-05  8.919518e-05   
MMT00082592   1.396383e-06  7.051967e-02  1.180095e-03  1.703197e-03   
MMT00082622   5.072608e-08  1.446964e-05  5.705016e-03  1.444670e-01   
MMT00082651   9.436231e-06  3.945290e-03  5.938962e-03  7.369466e-03   
MMT00082663   1.931199e-03  1.135110e-01  5.740455e-07  4.892323e-02   
MMT00082677   9.835944e-05  1.915738e-07  5.240328e-03  2.089488e-05   
MMT00082712   1.384538e-01  1.257881e-03  2.101181e-01  1.253344e-04   
MMT00082753   2.234119e-08  5.115783e-05  1.203741e-13  5.574431e-05   
MMT00082759   1.370091e-03  2.215884e-04  7.991345e-06  6.695043e-09   
MMT00082798   1.315196e-06  1.149304e+00  2.788908e-02  5.653087e-03   
MMT00082822   1.123074e-05  1.407546e-02  9.438444e-05  2.004508e-08   
MMT00082829   3.227641e-03  6.135732e-04  4.909215e-06  2.129642e-03   
MMT00082832   1.000000e+00  1.337377e-03  2.106797e-05  6.122326e-06   
MMT00082850   1.337377e-03  1.000000e+00  1.133054e-05  4.971757e-05   
MMT00082869   2.106797e-05  1.133054e-05  1.000000e+00  6.319144e-04   
MMT00082877   6.122326e-06  4.971757e-05  6.319144e-04  1.000000e+00   
MMT00082906   9.966405e-07  3.084718e-16  6.104880e-04  3.218101e-05   

substanceBXH   MMT00082906  
substanceBXH                
MMT00000044   3.406891e-06  
MMT00000046   1.795923e-07  
MMT00000051   1.692309e-18  
MMT00000080   1.166613e-01  
MMT00000102   3.554065e-12  
MMT00000149   8.391253e-04  
MMT00000159   2.011392e-07  
MMT00000207   1.362075e-01  
MMT00000212   2.771614e-07  
MMT00000241   3.077123e-03  
MMT00000268   2.745290e-03  
MMT00000283   1.483575e-05  
MMT00000334   5.419158e-12  
MMT00000373   8.178102e-12  
MMT00000401   1.404227e-11  
MMT00000418   6.578251e-03  
MMT00000464   1.157856e-11  
MMT00000517   1.617578e-05  
MMT00000549   3.038661e-06  
MMT00000602   1.859218e-05  
MMT00000608   1.449736e-07  
MMT00000701   2.252689e-04  
MMT00000713   4.029428e-06  
MMT00000743   1.479261e-05  
MMT00000793   1.831626e-09  
MMT00000840   2.564725e-07  
MMT00000887   5.186746e-07  
MMT00000963   1.618499e-13  
MMT00000988   6.457322e-11  
MMT00000996   1.510286e-06  
...                    ...  
MMT00082101   1.791036e-07  
MMT00082110   1.467396e-05  
MMT00082126   1.049279e-05  
MMT00082164   1.367593e-04  
MMT00082181   4.447977e-07  
MMT00082259   2.057165e-05  
MMT00082420   5.752435e-07  
MMT00082428   7.166605e-11  
MMT00082445   6.148024e-14  
MMT00082461   1.578867e-10  
MMT00082551   7.919985e-07  
MMT00082577   3.237007e-03  
MMT00082579   1.129107e-07  
MMT00082585   1.567330e-05  
MMT00082592   6.238053e-07  
MMT00082622   4.962094e-03  
MMT00082651   1.607230e-10  
MMT00082663   1.642063e-09  
MMT00082677   6.268290e-05  
MMT00082712   2.397598e-06  
MMT00082753   2.874204e-04  
MMT00082759   2.221897e-06  
MMT00082798   1.069492e-06  
MMT00082822   1.583379e-05  
MMT00082829   1.462919e-01  
MMT00082832   9.966405e-07  
MMT00082850   3.084718e-16  
MMT00082869   6.104880e-04  
MMT00082877   3.218101e-05  
MMT00082906   1.000000e+00  

[2785 rows x 2785 columns]```

```pythontest=np.float_power(temp,6)
np.fill_diagonal(test.values, 1.0)```

```pythonimport pandas as pd
import seaborn as sns  #用于绘制热图的工具包
from scipy.cluster import hierarchy  #用于进行层次聚类，话层次聚类图的工具包
from scipy import cluster   
import matplotlib.pyplot as plt
from sklearn import decomposition as skldec #用于主成分分析降维的包
```

```pythona=sns.clustermap(1-test.iloc[0:100,0:100],col_cluster=True,row_cluster=True,figsize=(25,25),cmap="seismic")```
*Output:*
```<Figure size 1800x1800 with 4 Axes>```

## 3、可视化模块（动态剪切树）

```pythonfrom dynamicTreeCut import cutreeHybrid
from scipy.spatial.distance import pdist
import numpy as np
from scipy.cluster.hierarchy import linkage,dendrogram

#d = np.transpose(np.arange(1,10001).reshape(100,100))
distances = pdist(1-test, "euclidean")```

```pythongeneTree=linkage(distances, "ward")
fig = plt.figure(figsize=(25, 10))
dn = dendrogram(geneTree)
dynamicMods=cutreeHybrid(geneTree,distM=distances,minClusterSize = 30,deepSplit = 2, pamRespectsDendro = False)```
*Output:*
```..cutHeight not given, setting it to 448.9691031625521  ===>  99% of the (truncated) height range in dendro.
```
```D:\Anaconda\lib\site-packages\pandas\core\series.py:842: FutureWarning: 
Passing list-likes to .loc or [] with any missing label will raise
KeyError in the future, you can use .reindex() as an alternative.

See the documentation here:
https://pandas.pydata.org/pandas-docs/stable/indexing.html#deprecate-loc-reindex-listlike
  return self.loc[key]
```
```..done.
```
```<Figure size 1800x720 with 1 Axes>```

```pythonset(dynamicMods['labels'])```
*Output:*
```{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}```

```pythonx=dn['ivl']
y=[dynamicMods['labels'][int(x)] for x in dn['ivl']]
yy=np.array([y])
```

```pythonz=[temp.index[int(x)] for x in dn['ivl']]
mol=pd.DataFrame(columns=['ivl','module','name'])
mol['ivl']=x
mol['module']=y
mol['name']=z
```

```pythonplt.pcolor(yy,cmap='seismic')```
*Output:*
```<matplotlib.collections.PolyCollection at 0x277a5cbd390>```
```<Figure size 432x288 with 1 Axes>```

```pythonfrom scipy.cluster import hierarchy
plt.figure(figsize=(25, 10))
grid = plt.GridSpec(3, 1, wspace=0.5, hspace=0.5)
#fig, (ax0, ax1) = plt.subplots(2, 1)
plt.subplot(grid[0:2,0])
hierarchy.set_link_color_palette(['#000000'])
ax0=hierarchy.dendrogram(geneTree,color_threshold=0, above_threshold_color='black')
plt.tick_params( \
    axis='x',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off')
#ax0.set_title('default: no edges')
plt.subplot(grid[2,0])
ax1=plt.pcolor(yy,cmap="seismic")
#ax1.set_title('thick edges')
#plt.subplots_adjust(hspace=0)
#fig.tight_layout()
plt.show()```
*Output:*
```D:\Anaconda\lib\site-packages\matplotlib\cbook\deprecation.py:107: MatplotlibDeprecationWarning: Passing one of 'on', 'true', 'off', 'false' as a boolean is deprecated; use an actual boolean (True/False) instead.
  warnings.warn(message, mplDeprecation, stacklevel=1)
D:\Anaconda\lib\site-packages\matplotlib\cbook\deprecation.py:107: MatplotlibDeprecationWarning: Passing one of 'on', 'true', 'off', 'false' as a boolean is deprecated; use an actual boolean (True/False) instead.
  warnings.warn(message, mplDeprecation, stacklevel=1)
D:\Anaconda\lib\site-packages\matplotlib\cbook\deprecation.py:107: MatplotlibDeprecationWarning: Passing one of 'on', 'true', 'off', 'false' as a boolean is deprecated; use an actual boolean (True/False) instead.
  warnings.warn(message, mplDeprecation, stacklevel=1)
```
```<Figure size 1800x720 with 2 Axes>```

```pythontemp.index```
*Output:*
```Index(['MMT00000044', 'MMT00000046', 'MMT00000051', 'MMT00000080',
       'MMT00000102', 'MMT00000149', 'MMT00000159', 'MMT00000207',
       'MMT00000212', 'MMT00000241',
       ...
       'MMT00082753', 'MMT00082759', 'MMT00082798', 'MMT00082822',
       'MMT00082829', 'MMT00082832', 'MMT00082850', 'MMT00082869',
       'MMT00082877', 'MMT00082906'],
      dtype='object', name='substanceBXH', length=2785)```

```pythondef all_list(arr):
    result = {}
    for i in set(arr):
        result[i] = arr.count(i)
    return result

# 结果：{0: 1, 1: 2, 2: 3, 3: 2}
testt=all_list(list(dynamicMods['labels']))```

```pythoni```
*Output:*
```'MMT00041434'```

```pythonann[ann['substanceBXH']=='MMT00041434']['gene_symbol'].iloc[0]```
*Output:*
```nan```

```python#增加基因名
ann=pd.read_csv('GeneAnnotation.csv')
rei=[ann[ann['substanceBXH']==i]['gene_symbol'].iloc[0] for i in mol['name']]
mol['gene']=rei```

```pythonmol[mol['module']==11].dropna().to_csv('mol11.csv')```

```pythonlen(mol)```
*Output:*
```2785```

## 5、主成分分析

```pythonpcamol=pd.DataFrame(columns=data.columns)
set_index=set(mol['module'])
for j in set_index:
    newdata=pd.DataFrame(columns=data.columns)
    for i in list(mol[mol['module']==j].dropna()['name']):
        newdata=newdata.append(data[data.index==i])
    from sklearn.decomposition import PCA
    pca = PCA(n_components=1) 
    reduced_X = pca.fit_transform(newdata.T)
    tepcamol=pd.DataFrame(reduced_X.T,columns=data.columns)
    pcamol=pcamol.append(tepcamol,ignore_index=True)
pcamol.index=set_index```

```python#性状矩阵
cri=pd.read_csv('ClinicalTraits.csv')
newcri=pd.DataFrame(columns=cri.columns)
for i in data.columns:   
    newcri=newcri.append(cri[cri['Mice']==i],ignore_index=True)
newcri```
*Output:*
```    Unnamed: 0    Mice Number Mouse_ID           Strain sex         DOB  \
0          207    F2_2      2    133-2  BxH ApoE-/-, F2   2  2001-10-22   
1          208    F2_3      3    133-3  BxH ApoE-/-, F2   2  2001-10-22   
2          219   F2_14     14    137-1  BxH ApoE-/-, F2   2  2001-10-22   
3          220   F2_15     15    137-2  BxH ApoE-/-, F2   2  2001-10-22   
4          224   F2_19     19    144-1  BxH ApoE-/-, F2   2  2001-10-24   
5          225   F2_20     20    144-2  BxH ApoE-/-, F2   2  2001-10-24   
6          228   F2_23     23    151-1  BxH ApoE-/-, F2   2  2001-11-10   
7          229   F2_24     24    151-2  BxH ApoE-/-, F2   2  2001-11-10   
8          231   F2_26     26    151-4  BxH ApoE-/-, F2   2  2001-11-10   
9          266   F2_37     37    155-2  BxH ApoE-/-, F2   2  2001-11-13   
10          87   F2_42     42    157-1  BxH ApoE-/-, F2   2  2001-11-15   
11          88   F2_43     43    157-2  BxH ApoE-/-, F2   2  2001-11-15   
12          90   F2_45     45    158-1  BxH ApoE-/-, F2   2  2001-11-15   
13          91   F2_46     46    158-2  BxH ApoE-/-, F2   2  2001-11-15   
14          92   F2_47     47    158-3  BxH ApoE-/-, F2   2  2001-11-15   
15          93   F2_48     48    158-4  BxH ApoE-/-, F2   2  2001-11-15   
16          96   F2_51     51    173-1  BxH ApoE-/-, F2   2  2001-12-04   
17          97   F2_52     52    173-2  BxH ApoE-/-, F2   2  2001-12-04   
18          99   F2_54     54    173-4  BxH ApoE-/-, F2   2  2001-12-04   
19         245   F2_63     63    178-1  BxH ApoE-/-, F2   2  2001-12-07   
20         247   F2_65     65    178-3  BxH ApoE-/-, F2   2  2001-12-07   
21         248   F2_66     66    178-4  BxH ApoE-/-, F2   2  2001-12-07   
22         250   F2_68     68    180-1  BxH ApoE-/-, F2   2  2001-12-10   
23         251   F2_69     69    180-2  BxH ApoE-/-, F2   2  2001-12-10   
24         252   F2_70     70    180-3  BxH ApoE-/-, F2   2  2001-12-10   
25         253   F2_71     71    181-1  BxH ApoE-/-, F2   2  2001-12-10   
26         254   F2_72     72    181-2  BxH ApoE-/-, F2   2  2001-12-10   
27         260   F2_78     78    198-2  BxH ApoE-/-, F2   2  2001-12-28   
28         261   F2_79     79    198-3  BxH ApoE-/-, F2   2  2001-12-28   
29         262   F2_80     80    198-4  BxH ApoE-/-, F2   2  2001-12-28   
..         ...     ...    ...      ...              ...  ..         ...   
105          1  F2_290    290    306-4  BxH ApoE-/-, F2   2  2002-03-22   
106          2  F2_291    291    307-1  BxH ApoE-/-, F2   2  2002-03-22   
107          7  F2_296    296    308-2  BxH ApoE-/-, F2   2  2002-03-22   
108          9  F2_298    298    308-4  BxH ApoE-/-, F2   2  2002-03-22   
109         10  F2_299    299    338-1  BxH ApoE-/-, F2   2  2002-03-22   
110         11  F2_300    300    338-2  BxH ApoE-/-, F2   2  2002-03-22   
111         13  F2_302    302    338-4  BxH ApoE-/-, F2   2  2002-03-22   
112         14  F2_303    303    373-1  BxH ApoE-/-, F2   2  2002-04-11   
113         15  F2_304    304    373-2  BxH ApoE-/-, F2   2  2002-04-11   
114         16  F2_305    305    373-3  BxH ApoE-/-, F2   2  2002-04-11   
115         17  F2_306    306    374-1  BxH ApoE-/-, F2   2  2002-04-11   
116         18  F2_307    307    374-2  BxH ApoE-/-, F2   2  2002-04-11   
117         19  F2_308    308    374-3  BxH ApoE-/-, F2   1  2002-04-14   
118         20  F2_309    309    327-1  BxH ApoE-/-, F2   2  2002-04-11   
119         21  F2_310    310    327-2  BxH ApoE-/-, F2   2  2002-04-11   
120         22  F2_311    311    327-3  BxH ApoE-/-, F2   2  2002-04-11   
121         23  F2_312    312    327-4  BxH ApoE-/-, F2   2  2002-04-11   
122        334  F2_320    320    343-1  BxH ApoE-/-, F2   2  2002-04-16   
123        335  F2_321    321    343-2  BxH ApoE-/-, F2   2  2002-04-16   
124        337  F2_323    323    343-4  BxH ApoE-/-, F2   2  2002-04-16   
125        338  F2_324    324    354-1  BxH ApoE-/-, F2   2  2002-05-06   
126        339  F2_325    325    354-2  BxH ApoE-/-, F2   2  2002-05-06   
127        340  F2_326    326    354-3  BxH ApoE-/-, F2   2  2002-05-06   
128        341  F2_327    327    354-4  BxH ApoE-/-, F2   2  2002-05-06   
129        342  F2_328    328    354-5  BxH ApoE-/-, F2   2  2002-05-06   
130        343  F2_329    329    359-1  BxH ApoE-/-, F2   2  2002-05-06   
131        344  F2_330    330    359-2  BxH ApoE-/-, F2   2  2002-05-06   
132        346  F2_332    332    359-4  BxH ApoE-/-, F2   2  2002-05-06   
133        186  F2_355    355    364-1  BxH ApoE-/-, F2   2  2002-05-31   
134        188  F2_357    357    364-3  BxH ApoE-/-, F2   2  2002-05-31   

      parents Western_Diet    Sac_Date       ...         Adiponectin  \
0       111.0   2001-12-17  2002-04-08       ...                 NaN   
1       111.0   2001-12-17  2002-04-08       ...              14.339   
2       109.0   2001-12-17  2002-04-08       ...              15.439   
3       109.0   2001-12-17  2002-04-08       ...              11.124   
4       110.0   2001-12-17  2002-04-08       ...              16.842   
5       110.0   2001-12-17  2002-04-08       ...              13.498   
6       112.0   2002-01-08  2002-04-30       ...              14.511   
7       112.0   2002-01-08  2002-04-30       ...              13.813   
8       112.0   2002-01-08  2002-04-30       ...              14.118   
9         NaN   2002-01-08  2002-04-30       ...              12.470   
10      111.0   2002-01-08  2002-04-30       ...              14.531   
11      111.0   2002-01-08  2002-04-30       ...               9.735   
12      111.0   2002-01-08  2002-04-30       ...               7.939   
13      111.0   2002-01-08  2002-04-30       ...              12.404   
14      111.0   2002-01-08  2002-04-30       ...              22.776   
15      111.0   2002-01-08  2002-04-30       ...              14.310   
16      109.0   2002-02-01  2002-05-17       ...              14.014   
17      109.0   2002-02-01  2002-05-17       ...              10.140   
18      109.0   2002-02-01  2002-05-17       ...               9.115   
19      112.0   2002-02-01  2002-05-17       ...               7.065   
20      112.0   2002-02-01  2002-05-17       ...               8.992   
21      112.0   2002-02-01  2002-05-17       ...              10.305   
22      111.0   2002-02-01  2002-05-17       ...               9.374   
23      111.0   2002-02-01  2002-05-17       ...               6.850   
24      111.0   2002-02-01  2002-05-17       ...               8.827   
25      111.0   2002-02-01  2002-05-17       ...                 NaN   
26      111.0   2002-02-01  2002-05-17       ...               9.239   
27      110.0   2002-02-26  2002-06-19       ...               5.767   
28      110.0   2002-02-26  2002-06-19       ...               4.676   
29      110.0   2002-02-26  2002-06-19       ...               4.820   
..        ...          ...         ...       ...                 ...   
105  229232.0   2002-05-14  2002-09-11       ...              11.274   
106     232.0   2002-05-14  2002-09-11       ...               7.099   
107     232.0   2002-05-14  2002-09-11       ...               8.710   
108     232.0   2002-05-14  2002-09-11       ...              14.237   
109     109.0   2002-05-14  2002-09-11       ...               7.592   
110     109.0   2002-05-14  2002-09-11       ...              14.366   
111     109.0   2002-05-14  2002-09-11       ...               8.404   
112       NaN   2002-07-24  2002-09-11       ...              12.048   
113       NaN   2002-07-24  2002-09-11       ...                 NaN   
114       NaN   2002-07-24  2002-09-11       ...              10.319   
115       NaN   2002-07-24  2002-09-11       ...              14.718   
116       NaN   2002-07-24  2002-09-11       ...              22.469   
117       NaN   2002-07-24  2002-09-11       ...              17.142   
118     230.0   2002-07-24  2002-11-13       ...              14.343   
119     230.0   2002-07-24  2002-11-13       ...              24.737   
120     230.0   2002-07-24  2002-11-13       ...              16.383   
121     230.0   2002-07-24  2002-11-13       ...              12.425   
122     228.0   2002-07-24  2002-11-13       ...              15.384   
123     228.0   2002-07-24  2002-11-13       ...              12.688   
124     228.0   2002-07-24  2002-11-13       ...               7.151   
125     228.0   2002-07-24  2002-11-13       ...               5.943   
126     228.0   2002-07-24  2002-11-13       ...              11.174   
127     228.0   2002-07-24  2002-11-13       ...              17.993   
128     228.0   2002-07-24  2002-11-13       ...              17.458   
129     228.0   2002-07-24  2002-11-13       ...               9.939   
130       NaN   2002-07-24  2002-11-13       ...              15.756   
131       NaN   2002-07-24  2002-11-13       ...               5.507   
132       NaN   2002-07-24  2002-11-13       ...               9.119   
133     228.0   2002-07-24  2002-11-13       ...              12.134   
134     228.0   2002-07-24  2002-11-13       ...                 NaN   

     Aortic lesions  Note  Aneurysm  Aortic_cal_M Aortic_cal_L  \
0          224500.0   NaN      56.0           5.0          0.0   
1          296250.0   NaN       8.0           4.0          NaN   
2          486313.0  0.75      27.0          12.0          NaN   
3          180750.0   NaN       0.0           0.0          NaN   
4          113000.0   NaN       0.0           0.0          NaN   
5          166750.0   NaN       6.0           0.0          NaN   
6          234000.0   NaN      28.0           8.0          NaN   
7          267500.0   NaN      33.0           8.0          NaN   
8          198000.0   NaN       0.0           0.0          0.0   
9          121000.0   NaN       0.0           0.0          0.0   
10         110000.0   NaN      42.0          14.0          0.0   
11         327250.0   NaN       0.0           0.0          0.0   
12          48250.0   NaN      18.0           1.0          0.0   
13         322000.0   NaN      25.0           1.0         10.0   
14         328250.0   NaN      42.0           3.0          0.0   
15         292500.0   NaN      56.0          10.0          0.0   
16         373250.0   NaN      34.0           2.0          5.0   
17         249000.0   NaN       0.0           0.0          0.0   
18          90500.0   NaN      24.0           6.0          0.0   
19         273250.0   NaN      48.0           4.0          0.0   
20         235250.0   NaN      33.0           0.0          0.0   
21         114500.0   NaN      23.0           6.0          0.0   
22         274750.0   NaN      39.0          15.0         21.0   
23         109250.0   NaN      12.0           2.0          3.0   
24         470500.0   NaN       4.0           0.0         17.0   
25         271750.0   NaN       5.0           4.0          3.0   
26         422250.0   NaN       3.0           0.0          0.0   
27          86000.0   NaN       6.0           4.0          0.0   
28         316500.0   NaN       0.0           2.0         11.0   
29         197500.0   NaN      18.0           0.0          0.0   
..              ...   ...       ...           ...          ...   
105        496250.0   NaN      16.0           0.0         17.0   
106             NaN   NaN      16.0           4.0          0.0   
107        133250.0   NaN       0.0           0.0          7.0   
108        301500.0   NaN       6.0           0.0          0.0   
109         57500.0   NaN       0.0           0.0          0.0   
110         62750.0   NaN       3.0           3.0          0.0   
111        126750.0   NaN      22.0          10.0          4.0   
112        154500.0   NaN      39.0           3.0          0.0   
113        352750.0   NaN       6.0           1.0          2.0   
114        190250.0   NaN      16.0           0.0          0.0   
115        166750.0   NaN       6.0           0.0          4.0   
116        275500.0   NaN       NaN           NaN          NaN   
117        192750.0   NaN      27.0           6.0          3.0   
118        318000.0   NaN      20.0           0.0          0.0   
119        268000.0   NaN      22.0           0.0          0.0   
120        183000.0   NaN       0.0           0.0          0.0   
121        411250.0   NaN       0.0           0.0         13.0   
122        238750.0   NaN      14.0           0.0          0.0   
123        312750.0   NaN      25.0           1.0         12.0   
124         83000.0   NaN       0.0           0.0          0.0   
125        129500.0   NaN      13.0           0.0          0.0   
126         57500.0   NaN      19.0           2.0          0.0   
127        246000.0   NaN      27.0          10.0          0.0   
128        313750.0   NaN      21.0           6.0          0.0   
129        221250.0   NaN       0.0           0.0          0.0   
130        143250.0   NaN       0.0           0.0         13.0   
131         60750.0   NaN       0.0           0.0          4.0   
132        366250.0   NaN      38.0           3.0          0.0   
133        109500.0   NaN      22.0           0.0          0.0   
134        214750.0   NaN      27.0           5.0          0.0   

     CoronaryArtery_Cal  Myocardial_cal  BMD_all_limbs  BMD_femurs_only  
0                   0.0             0.0            NaN              NaN  
1                   0.0             0.0            NaN              NaN  
2                   1.0             8.0            NaN              NaN  
3                   0.0             4.0            NaN              NaN  
4                   0.0             0.0            NaN              NaN  
5                   0.0             0.0            NaN              NaN  
6                   0.0             0.0            NaN              NaN  
7                   0.0             1.0            NaN              NaN  
8                   0.0             0.0            NaN              NaN  
9                   0.0             3.0            NaN              NaN  
10                  0.0             0.0            NaN              NaN  
11                  0.0             0.0            NaN              NaN  
12                  0.0             1.0            NaN              NaN  
13                  0.0             0.0            NaN              NaN  
14                  0.0             0.0            NaN              NaN  
15                  0.0             0.0            NaN              NaN  
16                  0.0             1.0            NaN              NaN  
17                  0.0             0.0            NaN              NaN  
18                  0.0             0.0            NaN              NaN  
19                  NaN             NaN            NaN              NaN  
20                  NaN             NaN            NaN              NaN  
21                  NaN             NaN            NaN              NaN  
22                  NaN             NaN            NaN              NaN  
23                  NaN             NaN            NaN              NaN  
24                  NaN             NaN            NaN              NaN  
25                  NaN             NaN            NaN              NaN  
26                  NaN             NaN            NaN              NaN  
27                  NaN             NaN            NaN              NaN  
28                  NaN             NaN            NaN              NaN  
29                  NaN             NaN            NaN              NaN  
..                  ...             ...            ...              ...  
105                 0.0             0.0            NaN              NaN  
106                 2.0             4.0         0.0548          0.07730  
107                 0.0             0.0         0.0621          0.09105  
108                 0.0             0.0            NaN              NaN  
109                 6.0             0.0         0.0575          0.08275  
110                 0.0            78.0         0.0523          0.07460  
111                 1.0            52.0         0.0588          0.08605  
112                 0.0             0.0            NaN              NaN  
113                 0.0             0.0         0.0553          0.07840  
114                 0.0             0.0            NaN              NaN  
115                 0.0             0.0            NaN              NaN  
116                 0.0             1.0            NaN              NaN  
117                 0.0             0.0         0.0502          0.07210  
118                 0.0             0.0         0.0573          0.07945  
119                 0.0             0.0         0.0560          0.08130  
120                 0.0             0.0         0.0565          0.07995  
121                 0.0             0.0         0.0551          0.07600  
122                 0.0             0.0         0.0591          0.08070  
123                 0.0             0.0         0.0556          0.08035  
124                 0.0             0.0         0.0588          0.08885  
125                 0.0             0.0         0.0529          0.07040  
126                 0.0             1.0         0.0520          0.07350  
127                 0.0             0.0         0.0525          0.07425  
128                 0.0             0.0         0.0509          0.06990  
129                 0.0             0.0         0.0529          0.07245  
130                 0.0             0.0         0.0507          0.06955  
131                 0.0             1.0         0.0588          0.08395  
132                 NaN             NaN         0.0582          0.08310  
133                 NaN             NaN         0.0564          0.07865  
134                 NaN             NaN            NaN              NaN  

[135 rows x 38 columns]```

```pythonnewnewcri=pd.DataFrame(columns=newcri.columns[10:14])
for i in range(len(newcri)):
    newnewcri=newnewcri.append(newcri.iloc[i,10:14])
newnewcri.index=newcri['Mice']
newnewcri.to_csv('newnewcri.csv')```

```pythonfrom scipy.stats import spearmanr,pearsonr,kendalltau
# seed random number generator
# calculate spearman's correlation
result_1=pd.DataFrame(columns=newnewcri.columns)
result_p=pd.DataFrame(columns=newnewcri.columns)
for j in newnewcri.columns:
    co=[]
    pvv=[]
    for i in range(len(pcamol)):   
        tempcor=pd.DataFrame(columns=['x','y'])
        tempcor['x']=list(newcri[j])
        tempcor['y']=list(pcamol.iloc[i])
        tempcor=tempcor.dropna()
        coef,pv=pearsonr(tempcor['x'],tempcor['y'])
        co.append(coef)
        pvv.append(pv)
    result_1[j]=co
    result_p[j]=pvv
        #print(coef)
result_1=abs(result_1)```

```pythonstr(coef)+'('+str(pv)+')'```
*Output:*
```'-0.4465721458538156(5.667571913912509e-08)'```

```pythonresult_1.index=set_index
result_p.index=set_index
result_p```
*Output:*
```        weight_g  length_cm        ab_fat     other_fat
1   9.833129e-01   0.123533  7.370441e-01  1.315513e-02
2   1.314628e-05   0.339842  5.819885e-04  3.043035e-05
3   1.126663e-03   0.805638  1.221870e-03  1.409295e-02
4   8.569493e-02   0.161288  1.807469e-02  5.889726e-01
5   2.259025e-03   0.068516  3.314510e-04  5.583534e-01
6   9.331135e-01   0.428312  7.461927e-01  9.601241e-01
7   9.833476e-01   0.501911  7.355670e-01  9.706715e-02
8   5.565951e-04   0.113775  1.661000e-03  1.228199e-01
9   1.512899e-02   0.109674  2.896697e-03  9.617716e-01
10  5.794533e-02   0.705549  1.789780e-01  4.638594e-02
11  9.573529e-19   0.040942  9.197359e-12  5.490529e-11
12  4.451714e-13   0.198976  1.392220e-09  5.667572e-08```

```pythonplt.figure(figsize=(10,10))
plt.rc('font', family='Arial')
plt.rcParams["font.weight"] = "bold"
plt.rcParams["font.size"] = "12"
plt.rcParams["axes.labelweight"] = "bold"
sns.heatmap(result_1,vmin=0, vmax=1,cmap='YlGnBu',annot=True,square=True)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)```
*Output:*
```(array([ 0.5,  1.5,  2.5,  3.5,  4.5,  5.5,  6.5,  7.5,  8.5,  9.5, 10.5,
        11.5]), <a list of 12 Text yticklabel objects>)```
```<Figure size 720x720 with 2 Axes>```

```python#!/usr/bin/env python

##############################################################
## The script prints out the p-value of STRING protein-protein
## interaction enrichment method for the given set of proteins 
##
## Requires requests module:
## type "python -m pip install requests" in command line (win)
## or terminal (mac/linux) to install the module
##############################################################

import requests ## python -m pip install requests

string_api_url = "https://string-db.org/api"
output_format = "tsv-no-header"
method = "ppi_enrichment"

##
## Construct the request
##

request_url = "/".join([string_api_url, output_format, method])

##
## Set parameters
##

my_genes = ['7227.FBpp0074373', '7227.FBpp0077451', '7227.FBpp0077788',
            '7227.FBpp0078993', '7227.FBpp0079060', '7227.FBpp0079448']

params = {

    "identifiers" : "%0d".join(my_genes), # your proteins
    "species" : 7227, # species NCBI identifier 
    "caller_identity" : "www.awesome_app.org" # your app name

}

##
## Call STRING
##

response = requests.post(request_url, data=params)

##
## Parse and print the respons Parse and print the responsee
##

for line in response.text.strip().split("\n"):
    pvalue = line.split("\t")[5]
    print("P-value:", pvalue)```
*Output:*
```P-value: 6.34e-12
```

```pythonmol[mol['module']==11].dropna()['gene'].values```
*Output:*
```array(['2210415F13Rik', 'AI747448', 'D630035O19Rik', ..., 'Stat1', 'Ccr5',
       'Skil'], dtype=object)```

```pythonimport gseapy as gp```
*Output:*
```Creating directory C:\Users\FernandoZeng\AppData\Local\bioservices\bioservices 
```

```pythonmol[mol['module']==11].dropna()['gene'].to_csv('gene_list.txt',index=False)```

```pythonenr = gp.enrichr(gene_list=mol[mol['module']==11].dropna()['gene'].values.tolist(),
                 gene_sets=['KEGG_2016'],
                 organism='Human', # don't forget to set organism to the one you desired! e.g. Yeast
                 description='test_name',
                 outdir='test/enrichr_kegg',
                 # no_plot=True,
                 cutoff=0.5 # test dataset, use lower value from range(0,1)
                )```

```pythonenr.results.head(5)```
*Output:*
```    Gene_set                                               Term Overlap  \
0  KEGG_2016  AGE-RAGE signaling pathway in diabetic complic...  18/101   
1  KEGG_2016  Cytokine-cytokine receptor interaction Homo sa...  34/265   
2  KEGG_2016            Carbon metabolism Homo sapiens hsa01200  19/113   
3  KEGG_2016  Biosynthesis of amino acids Homo sapiens hsa01230   14/74   
4  KEGG_2016  Amino sugar and nucleotide sugar metabolism Ho...   10/48   

    P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  Odds Ratio  \
0  0.000051          0.006630            0                     0    3.274594   
1  0.000058          0.006630            0                     0    2.233731   
2  0.000073          0.006630            0                     0    3.052700   
3  0.000178          0.012031            0                     0    3.516183   
4  0.000664          0.035964            0                     0    3.957522   

   Combined Score                                              Genes  
0       32.342930  STAT5A;JUN;EDN1;VCAM1;PRKCB;STAT1;SERPINE1;STA...  
1       21.809205  CSF3R;FLT1;FLT4;PDGFA;IL1RAP;TNF;CSF2RA;EGFR;C...  
2       29.060649  SDS;GOT1;GPT2;IDH1;PGAM2;ENO1;PGD;ENO3;MUT;ACA...  
3       30.366033  SDS;GOT1;TAT;GPT2;IDH1;PGAM2;ENO1;ENO3;ASS1;CT...  
4       28.960780   UGDH;GALE;GALT;UGP2;HEXB;GNPNAT1;HK2;GCK;HK1;GNE  ```

```pythonenr.res2d```
*Output:*
```      Gene_set                                               Term Overlap  \
0    KEGG_2013                  HSA00100 BIOSYNTHESIS OF STEROIDS    9/24   
1    KEGG_2013    HSA04060 CYTOKINE CYTOKINE RECEPTOR INTERACTION  33/257   
2    KEGG_2013              HSA00520 NUCLEOTIDE SUGARS METABOLISM     4/6   
3    KEGG_2013          HSA00252 ALANINE AND ASPARTATE METABOLISM    8/33   
4    KEGG_2013                        HSA05010 ALZHEIMERS DISEASE    7/28   
5    KEGG_2013                            HSA04510 FOCAL ADHESION  24/200   
6    KEGG_2013                 HSA00521 STREPTOMYCIN BIOSYNTHESIS    4/10   
7    KEGG_2013                      HSA00052 GALACTOSE METABOLISM    7/32   
8    KEGG_2013                    HSA03320 PPAR SIGNALING PATHWAY   11/70   
9    KEGG_2013                    HSA00900 TERPENOID BIOSYNTHESIS     3/6   
10   KEGG_2013           HSA04920 ADIPOCYTOKINE SIGNALING PATHWAY   11/72   
11   KEGG_2013           HSA00330 ARGININE AND PROLINE METABOLISM    7/35   
12   KEGG_2013      HSA04620 TOLL LIKE RECEPTOR SIGNALING PATHWAY  13/102   
13   KEGG_2013                   HSA00401 NOVOBIOCIN BIOSYNTHESIS     2/3   
14   KEGG_2013    HSA00760 NICOTINATE AND NICOTINAMIDE METABOLISM    5/24   
15   KEGG_2013                HSA04630 JAK STAT SIGNALING PATHWAY  17/153   
16   KEGG_2013                 HSA04930 TYPE II DIABETES MELLITUS    7/44   
17   KEGG_2013                  HSA04614 RENIN ANGIOTENSIN SYSTEM    4/17   
18   KEGG_2013                  HSA04512 ECM RECEPTOR INTERACTION   11/87   
19   KEGG_2013   HSA00260 GLYCINE SERINE AND THREONINE METABOLISM    7/45   
20   KEGG_2013                        HSA00310 LYSINE DEGRADATION    7/47   
21   KEGG_2013                    HSA00530 AMINOSUGARS METABOLISM    5/29   
22   KEGG_2013                    HSA00480 GLUTATHIONE METABOLISM    6/39   
23   KEGG_2013                   HSA00950 ALKALOID BIOSYNTHESIS I     2/5   
24   KEGG_2013                         HSA00230 PURINE METABOLISM  15/145   
25   KEGG_2013                    HSA00565 ETHER LIPID METABOLISM    5/31   
26   KEGG_2013            HSA00010 GLYCOLYSIS AND GLUCONEOGENESIS    8/64   
27   KEGG_2013                   HSA00061 FATTY ACID BIOSYNTHESIS     2/6   
28   KEGG_2013   HSA01040 POLYUNSATURATED FATTY ACID BIOSYNTHESIS    3/14   
29   KEGG_2013                             HSA05060 PRION DISEASE    3/14   
..         ...                                                ...     ...   
130  KEGG_2013                       HSA00350 TYROSINE METABOLISM    3/59   
131  KEGG_2013     HSA04070 PHOSPHATIDYLINOSITOL SIGNALING SYSTEM    4/78   
132  KEGG_2013                      HSA00340 HISTIDINE METABOLISM    2/41   
133  KEGG_2013                                  HSA03010 RIBOSOME    5/98   
134  KEGG_2013                            HSA05219 BLADDER CANCER    2/42   
135  KEGG_2013          HSA00532 CHONDROITIN SULFATE BIOSYNTHESIS    1/22   
136  KEGG_2013                 HSA04020 CALCIUM SIGNALING PATHWAY   9/174   
137  KEGG_2013   HSA00361 GAMMA HEXACHLOROCYCLOHEXANE DEGRADATION    1/23   
138  KEGG_2013  HSA00563 GLYCOSYLPHOSPHATIDYLINOSITOL ANCHOR B...    1/23   
139  KEGG_2013                  HSA04940 TYPE I DIABETES MELLITUS    2/45   
140  KEGG_2013                             HSA04360 AXON GUIDANCE   6/128   
141  KEGG_2013                     HSA04310 WNT SIGNALING PATHWAY   7/149   
142  KEGG_2013  HSA04650 NATURAL KILLER CELL MEDIATED CYTOTOXI...   6/132   
143  KEGG_2013                        HSA05213 ENDOMETRIAL CANCER    2/52   
144  KEGG_2013           HSA00903 LIMONENE AND PINENE DEGRADATION    1/29   
145  KEGG_2013             HSA01032 GLYCAN STRUCTURES DEGRADATION    1/30   
146  KEGG_2013                   HSA04140 REGULATION OF AUTOPHAGY    1/30   
147  KEGG_2013                  HSA05220 CHRONIC MYELOID LEUKEMIA    3/76   
148  KEGG_2013                       HSA05040 HUNTINGTONS DISEASE    1/31   
149  KEGG_2013                      HSA05217 BASAL CELL CARCINOMA    2/56   
150  KEGG_2013  HSA04130 SNARE INTERACTIONS IN VESICULAR TRANS...    1/36   
151  KEGG_2013                 HSA00190 OXIDATIVE PHOSPHORYLATION   5/128   
152  KEGG_2013                HSA01510 NEURODEGENERATIVE DISEASES    1/38   
153  KEGG_2013            HSA04120 UBIQUITIN MEDIATED PROTEOLYSIS    1/40   
154  KEGG_2013      HSA00860 PORPHYRIN AND CHLOROPHYLL METABOLISM    1/41   
155  KEGG_2013          HSA01030 GLYCAN STRUCTURES BIOSYNTHESIS 1   4/117   
156  KEGG_2013             HSA00562 INOSITOL PHOSPHATE METABOLISM    1/51   
157  KEGG_2013                        HSA04742 TASTE TRANSDUCTION    1/53   
158  KEGG_2013       HSA04612 ANTIGEN PROCESSING AND PRESENTATION    2/83   
159  KEGG_2013   HSA04080 NEUROACTIVE LIGAND RECEPTOR INTERACTION   9/254   

      P-value  Adjusted P-value  Old P-value  Old Adjusted P-value  \
0    0.000008          0.001288            0                     0   
1    0.000072          0.005783            0                     0   
2    0.000208          0.011109            0                     0   
3    0.000791          0.031658            0                     0   
4    0.001382          0.044221            0                     0   
5    0.001708          0.045554            0                     0   
6    0.002379          0.054377            0                     0   
7    0.003154          0.063089            0                     0   
8    0.004015          0.068166            0                     0   
9    0.004260          0.068166            0                     0   
10   0.005004          0.071327            0                     0   
11   0.005350          0.071327            0                     0   
12   0.011245          0.129064            0                     0   
13   0.011293          0.129064            0                     0   
14   0.015032          0.155651            0                     0   
15   0.015565          0.155651            0                     0   
16   0.018714          0.167174            0                     0   
17   0.018969          0.167174            0                     0   
18   0.019852          0.167174            0                     0   
19   0.020999          0.167996            0                     0   
20   0.026137          0.199137            0                     0   
21   0.032467          0.230610            0                     0   
22   0.033169          0.230610            0                     0   
23   0.034591          0.230610            0                     0   
24   0.038640          0.247298            0                     0   
25   0.041939          0.258087            0                     0   
26   0.045873          0.271838            0                     0   
27   0.049757          0.283979            0                     0   
28   0.053246          0.283979            0                     0   
29   0.053246          0.283979            0                     0   
..        ...               ...          ...                   ...   
130  0.724155          0.883358            0                     0   
131  0.728771          0.883358            0                     0   
132  0.737159          0.886807            0                     0   
133  0.743330          0.887558            0                     0   
134  0.749257          0.888009            0                     0   
135  0.759566          0.891693            0                     0   
136  0.769619          0.891693            0                     0   
137  0.774658          0.891693            0                     0   
138  0.774658          0.891693            0                     0   
139  0.782711          0.894527            0                     0   
140  0.820791          0.929780            0                     0   
141  0.832331          0.929780            0                     0   
142  0.842322          0.929780            0                     0   
143  0.845896          0.929780            0                     0   
144  0.847281          0.929780            0                     0   
145  0.856870          0.929780            0                     0   
146  0.856870          0.929780            0                     0   
147  0.863156          0.929780            0                     0   
148  0.865858          0.929780            0                     0   
149  0.874018          0.932286            0                     0   
150  0.903011          0.956667            0                     0   
151  0.909740          0.956667            0                     0   
152  0.914813          0.956667            0                     0   
153  0.925180          0.959876            0                     0   
154  0.929880          0.959876            0                     0   
155  0.940508          0.964624            0                     0   
156  0.963358          0.975969            0                     0   
157  0.967820          0.975969            0                     0   
158  0.969870          0.975969            0                     0   
159  0.980705          0.980705            0                     0   

     Odds Ratio  Combined Score  \
0      9.026988      105.882816   
1      2.234797       21.308560   
2     29.990400      254.214946   
3      4.807961       34.336695   
4      5.005346       32.956693   
5      2.058758       13.118986   
6      9.994667       60.378520   
7      4.203593       24.208269   
8      2.802910       15.465954   
9     14.982414       81.779948   
10     2.710721       14.360041   
11     3.752606       19.628921   
12     2.195955        9.855127   
13    29.944089      134.256327   
14     3.945683       16.562169   
15     1.880558        7.828244   
16     2.838445       11.292677   
17     4.611200       18.283212   
18     2.173964        8.520749   
19     2.763601       10.676508   
20     2.625140        9.567101   
21     3.122832       10.703557   
22     2.726253        9.286032   
23     9.980298       33.575210   
24     1.733656        5.640379   
25     2.882306        9.141335   
26     2.142857        6.604045   
27     7.484824       22.458994   
28     4.084369       11.978761   
29     4.084369       11.978761   
..          ...             ...   
130    0.800360        0.258315   
131    0.807438        0.255470   
132    0.766241        0.233667   
133    0.802922        0.238159   
134    0.747045        0.215651   
135    0.711625        0.195703   
136    0.814064        0.213171   
137    0.679243        0.173434   
138    0.679243        0.173434   
139    0.694814        0.170223   
140    0.733922        0.144940   
141    0.735444        0.134972   
142    0.710470        0.121912   
143    0.597316        0.099966   
144    0.533520        0.088416   
145    0.515095        0.079566   
146    0.515095        0.079566   
147    0.613416        0.090270   
148    0.497898        0.071715   
149    0.552952        0.074458   
150    0.426656        0.043527   
151    0.606111        0.057336   
152    0.403550        0.035930   
153    0.382815        0.029770   
154    0.373224        0.027133   
155    0.527660        0.032364   
156    0.298420        0.011140   
157    0.286911        0.009385   
158    0.368102        0.011262   
159    0.545886        0.010636   

                                                 Genes  
0       FDPS;IDI1;SQLE;NSDHL;PMVK;MVD;DHCR7;LSS;SC4MOL  
1    CSF3R;FLT1;FLT4;IL1RAP;TNF;CSF2RA;EGFR;CXCL5;I...  
2                                  UGDH;GALE;GALT;UGP2  
3               GOT1;GPT2;ASL;DLAT;AGXT;ASPA;CRAT;ASS1  
4                   C1QB;C1QA;BACE2;MME;CASP3;IL1B;TNF  
5    JUN;FLT1;LAMB3;CAV2;LAMA1;TNC;PDGFA;PRKCA;PARV...  
6                                    IMPA2;GCK;HK2;HK1  
7                     GALE;GALT;UGP2;RDH11;HK2;GCK;HK1  
8    FADS2;FABP3;ACSL1;EHHADH;FABP7;DBI;ACSL4;PPARG...  
9                                       FDPS;IDI1;SQLE  
10   SOCS3;PRKAA2;ACSL1;FRAP1;STAT3;ACSL4;CD36;JAK2...  
11                   GOT1;P4HA1;ASL;CKB;PRODH;ASS1;OTC  
12   JUN;STAT1;PIK3R1;TNF;CXCL10;IL1B;CCL5;CCL4;IRF...  
13                                            GOT1;TAT  
14                         NT5E;NNMT;PBEF1;AOX1;NUDT12  
15   STAT5A;CSF3R;STAT1;STAT3;PIK3R1;OSMR;PRLR;CSF2...  
16               SOCS3;SOCS1;FRAP1;PIK3R1;MAFA;TNF;GCK  
17                                  ENPEP;MME;CTSG;AGT  
18   VTN;LAMB3;COL4A1;LAMA1;COL5A3;CHAD;TNC;ITGA8;S...  
19                 ALAS2;SDS;RDH11;CTH;PIPOX;AGXT;GNMT  
20            TMLHE;RDH11;NSD1;EHHADH;PIPOX;ACAT2;AASS  
21                            HEXB;GNPNAT1;HK2;HK1;GNE  
22                   GSTM2;GSTM1;MGST3;IDH1;GSTA2;GPX7  
23                                            GOT1;TAT  
24   GUCY1A3;PDE6H;AK1;NME2;ADCY1;FHIT;PAICS;DCK;PA...  
25                 ENPP2;PLA2G6;AGPAT2;AGPAT4;PAFAH1B1  
26               PGAM2;GALM;ENO1;DLAT;ENO3;HK2;GCK;HK1  
27                                          FASN;ACACB  
28                                    FADS2;FASN;FADS1  
29                                    LAMA1;TNF;NFE2L2  
..                                                 ...  
130                                      GOT1;TAT;AOX1  
131                           IMPA2;PRKCA;ITPR3;PIK3R1  
132                                           HAL;ASPA  
133                      RPL31;RPS8;RPL36A;RPL29;RPL3L  
134                                        RASSF1;EGFR  
135                                              CHSY1  
136  CYSLTR1;EDNRB;ERBB3;PRKCA;ITPR3;CACNA1S;AVPR1A...  
137                                               PON2  
138                                              GPLD1  
139                                           IL1B;TNF  
140                 SEMA5A;UNC5B;SEMA3D;FYN;NGEF;GNAI1  
141             FZD3;JUN;CCND2;PRKCA;CSNK1E;TBL1X;NKD1  
142                   CASP3;PRKCA;FYN;PIK3R1;TNF;CD244  
143                                        PIK3R1;EGFR  
144                                             EHHADH  
145                                               HEXB  
146                                             PRKAA2  
147                               STAT5A;PIK3R1;TGFBR2  
148                                              CASP3  
149                                          FZD3;BMP2  
150                                               BET1  
151                NDUFA5;COX7A2;ATP6V1E1;COX7A1;COX7C  
152                                              CASP3  
153                                               CUL2  
154                                              ALAS2  
155                         CHSY1;GALNT3;GCNT1;C1GALT1  
156                                              IMPA2  
157                                              ITPR3  
158                                           B2M;CTSB  
159  GABBR1;CYSLTR1;MTNR1A;EDNRB;HRH4;GLP2R;CTSG;AV...  

[160 rows x 10 columns]```

```python# simple plotting function
from gseapy.plot import barplot, dotplot

# to save your figure, make sure that ``ofname`` is not None
barplot(enr.res2d,title='KEGG_2016',)```
*Output:*
```<Figure size 468x432 with 1 Axes>```

```python# to save your figure, make sure that ``ofname`` is not None
dotplot(enr.res2d, title='KEGG_2013',cmap='viridis_r',cutoff=10,ofname='kegg.png')```

```python
module_list=mol[mol['module']==11].dropna()['gene'].values.tolist()
#enrichment of KEGG and GO
module_KEGG=enrichment_KEGG(module_list,
                     organism='Mouse')
module_GO=enrichment_GO(module_list,
                     organism='Mouse')
#enrichment of GSEA
result_DEG=find_DEG(data)
module_GSEA=enrichment_GSEA(data=result_DEG,
               gene_sets='KEGG_2016',
               processes=4,
               permutation_num=100)```

```pythonk=dotplot(ttty, title='go',cmap='viridis_r')
print(k)```
*Output:*
```AxesSubplot(0.125,0.125;0.775x0.755)
```
```<Figure size 432x396 with 2 Axes>```

```pythondef enrichment_GO(gene_list,
                    go_mode='Bio',
                    organism='Human',
                    description='test_name',
                    outdir='enrichment_go',
                    cutoff=0.5):
    if(go_mode=='Bio'):
        geneset='GO_Biological_Process_2018'
    if(go_mode=='Cell'):
        geneset='GO_Cellular_Component_2018'
    if(go_mode=='Mole'):
        geneset='GO_Molecular_Function_2018'
    enr = gp.enrichr(gene_list=gene_list,
                 gene_sets=geneset,
                 organism=organism, # don't forget to set organism to the one you desired! e.g. Yeast
                 description=description,
                 outdir=outdir,
                 # no_plot=True,
                 cutoff=cutoff # test dataset, use lower value from range(0,1)
                )
    dotplot(enr.res2d, title=description,cmap='viridis_r',ofname='go_'+go_mode+'.png')
    return enr.res2d```

```python#
result_DEG=find_DEG(data)
enrichment_GSEA(data=result_DEG,
               gene_sets='KEGG_2016',
               processes=4,
               permutation_num=100)```

```pythondef enrichment_GSEA(data,
                   gene_sets='KEGG_2016',
                   processes=4
                   permutation_num=100,
                   outdir='prerank_report_kegg',
                   seed=6):
    rnk=pd.DataFrame(columns=['genename','FoldChange'])
    rnk['genename']=data.index
    rnk['FoldChange']=data['FoldChange'].tolist()
    rnk1=rnk1.drop_duplicates(['genename'])
    rnk1=rnk1.sort_values(by='FoldChange', ascending=False)
    
    pre_res = gp.prerank(rnk=rnk1, gene_sets=gene_sets,
                     processes=processes,
                     permutation_num=permutation_num, # reduce number to speed up testing
                     outdir=outdir, format='png', seed=seed)
    pre_res.res2d.sort_index().to_csv('GSEA_result.csv')
    return pre_res```

```pythondef density_norm(data):
    grid = plt.GridSpec(1, 4, wspace=1, hspace=0.1)
    #plt.subplot(grid[0,0:2])
    data.plot(kind = 'density')
    plt.title('Pre')
    #plt.xlim(-1000,1000)
    ERlist=ERgene.FindERG(data)
    data2=ERgene.normalizationdata(data,ERlist[0])
    #plt.subplot(grid[0,2:4])
    data2.plot(kind='density')
    plt.title('After')```

