# Eular

"""
    read x0, y0, h, x_target 
    n <- round((x_target - x0)/h) 
    x <- x0; y <- y0 
    repeat n times: y <- y + h*f(x,y); x <- x + h 
    print y
"""
def eular_method(f,x0,y0,h,target) :
    if h <= 0:
        print("Error: Step size h must be greater than 0.")
        return None

    if target < x0:
        print("Error: Target must be greater than or equal to x0.")
        return None
    x = x0
    y = y0
    print(f"{"Steps":<6} | {"x":<10} | {"y":<10}")
    print("---------------------------------------")
    print(f"{"0":<6} | {x0:<10.2f} | {y0:<10.4f}")
    steps = int((target - x0) / h)
    for i in range(1,steps+1) :
        y = y + h * f(x,y)
        x = x + h
        print(f"{i:<6} | {x:<10.2f} | {y:<10.4f}")
    return y

if __name__ == "__main__" :
    f = lambda x,y : -2.2067e-12 * (y**4-81e8)
    x0 = 0
    y0 = 1200
    h = 60
    target_x = 480
    value_y = eular_method(f,x0,y0,h,target_x)

    print(f"\nFinal Estimated y at x = {target_x} is {value_y:.4f}")

"""
How many integration steps are required, and how is that number computed? 
Answer: n=(x_target-x0)/h=(1.0-0)/0.2=5 steps. 
2. Show the numerical quantities used in the first update and report the first computed y-value. 
Answer: At the first step, f(0,0.5)=1.500000; therefore y1=0.500000+0.2(1.500000)=0.800000. 
3. In a Python implementation, what loop condition or loop count prevents stepping beyond the target x-value? 
Answer: Use exactly n=5 iterations (for example, for _ in range(n)) after checking that (x_target-x0)/h is an integer within tolerance.
What is the global order of accuracy of Euler's method, and what does this imply when h is approximately halved? 
Answer: Euler's method is first-order globally. Halving h typically reduces the leading global error by about a factor of 2 for a 
sufficiently smooth problem. 
 State one programming check that should be included before running this method. 
Answer: Check h>0, x_target>=x0, and that the interval length is divisible by h within a small floating-point tolerance; also ensure f(x,y) 
returns a finite numeric value. 
"""

# RK-2 Method :

"""
    Heun's Method :
    read inputs; compute n; set x=x0, y=y0 
    repeat n times: 
    k1 <- f(x,y) 
    k2 <- f(x+h, y+h*k1) 
    y <- y + h*(k1+k2)/2; x <- x+h 
    print y 
    MidPoint :
    read inputs; compute n; set x=x0, y=y0 
    repeat n times: 
    k1 <- f(x,y) 
    k2 <- f(x+h/2, y+h*k1/2) 
    y <- y + h*k2; x <- x+h 
    print y
    Ralston's Method :
    read inputs; compute n; set x=x0, y=y0 
    repeat n times: 
    k1 <- f(x,y) 
    k2 <- f(x+2h/3, y+2h*k1/3) 
    y <- y + h*(k1/4 + 3*k2/4); x <- x+h 
    print y 
"""
def rk2_general(f,x0,y0,h,target,a1,a2,p1,q11) :
    # Error handling
    if h <= 0:
        print("Error: Step size h must be greater than 0.")
        return None

    if target < x0:
        print("Error: Target value must be greater than or equal to x0.")
        return None

    if a1 < 0 or a2 < 0:
        print("Error: a1 and a2 must be non-negative.")
        return None

    x = x0
    y = y0
    steps = int((target - x0) /h)
    print(f"{"steps":<6} | {"x":<10} | {"y":<10} | {"k1":<10} | {"k2":<10}")
    print("---------------------------------------------------------------")
    print(f"{0:<6} | {x:<10.2f} | {y:<10.4f} | {'........':<10} | {'........':<10}")
    for i in range(1,steps+1) :
        k1 = f(x,y)
        k2 = f(x + p1*h,y + q11*k1*h)
        y = y + (a1*k1 + a2*k2) * h
        x = x + h
        print(f"{i:<6} | {x:<10.2f} | {y:<10.4f} | {k1:<10.4f} | {k2:<10.4f}")
if __name__ == "__main__" :
    f = lambda x,y : -2.2067e-12 * (y**4 - 81e8)
    x0 = 0
    y0 = 1200
    h = 240
    target_x = 480 

    print("---------- Heun's Method (a2 = 1/2) ----------")
    rk2_general(f,x0,y0,h,target_x,a1=0.5,a2=0.5,p1=1.0,q11=1.0)
    print("\n============================================================")
    print("\n---------- Midpoint Method (a2 = 1) ----------")
    rk2_general(f,x0,y0,h,target_x,a1=0.0,a2=1.0,p1=0.5,q11=0.5)
    print("\n============================================================")
    print("\n----------Ralston Method (a2 = 2/3) ----------")
    rk2_general(f,x0,y0,h,target_x,a1=1/3,a2=2/3,p1=3/4,q11=3/4)

'''
Show the numerical quantities used in the first update and report the first computed y-value. 
Answer: k1=-1.000000, predictor=0.800000, k2=-0.601331; y1=0.839867.

What is the global order of accuracy of this RK2 method, and what does this imply when h is approximately halved? 
Answer: The method is second-order globally. Halving h typically reduces the leading global error by about a factor of 4 for a 
sufficiently smooth problem. 
'''

# RK-4 Mthod :

"""
read inputs; compute n; set x=x0, y=y0 
repeat n times: 
  k1 <- f(x,y) 
  k2 <- f(x+h/2, y+h*k1/2) 
  k3 <- f(x+h/2, y+h*k2/2) 
  k4 <- f(x+h, y+h*k3) 
  y <- y + h*(k1+2*k2+2*k3+k4)/6; x <- x+h 
print y 
"""

def r4_generate(f,x0,y0,h,target) :
    # Error handling
    if h <= 0:
        print("Error: Step size h must be greater than 0.")
        return None
    
    if target < x0:
        print("Error: Target value must be greater than or equal to x0.")
        return None
    
    x = x0
    y = y0
    steps = int((target - x0) / h)
    print(f"{"Steps":<6} | {"x":<6} | {"y":<6} | {"k1":<6} | {"k2":<6} | {"k3":<6} | {"k4":<6}")
    print("--------------------------------------------------------------------------------------")
    print(f"{0:<6} | {x:<6.2f} | {y:<6.4f} | {"......":<6} | {"......":<6} | {"......":<6} | {"......":<6}")
    for i in range(1,steps+1) :
        k1 = f(x,y)
        k2 = f(x + 0.5*h, y + 0.5*k1*h)
        k3 = f(x + 0.5*h, y + 0.5*k2*h)
        k4 = f(x + h, y + k3*h)

        y = y + (h * (k1 + 2*k2 + 2*k3 + k4)) / 6
        x = x + h
        print(f"{i:<6} | {x:<6.2f} | {y:<6.4f} | {k1:<6.4f} | {k2:<6.4f} | {k3:<6.4f} | {k4:<6.4f}")
    return y

if __name__ == "__main__" :
    f = lambda x, y: -2.2067e-12 * (y**4 - 81e8)
    x0 = 0
    y0 = 1200
    h = 240
    target_x = 480
    value_y = r4_generate(f,x0,y0,h,target_x)
    print(f"Final Estimated y at x = {target_x} is {value_y:.6f}")

"""
What is the global order of accuracy of classical RK4, and what does this imply when h is approximately halved? 
Answer: RK4 is fourth-order globally. Halving h typically reduces the leading global error by about a factor of 16 for a sufficiently 
smooth problem. 
5. State one programming check that should be included before running this method. 
Answer: Check h>0, x_target>=x0, and that the interval length is divisible by h within a small floating-point tolerance; also ensure f(x,y) 
returns a finite numeric value.
"""