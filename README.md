# Lazy matrix
lazy_matrix defines a set of template classes to operate on matrices lazily. Only the elements
explicitely required are computed.  
For instance, the following code only evaluates m1(1,1) + m2(1,1):

```
lazy_matrix<double> m1,m2;
m1.init(3,3);
m2.init(3,3);

auto m3 = m1 + m2

m1.setValue(1,1,2.0);
m2.setValue(1,1,12.0);

std::cout << m3(1,1) << std::endl;
```

The output will be 14.0.
