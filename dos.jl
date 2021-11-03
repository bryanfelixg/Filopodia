using BumpFuncs
using Plots
##

function two_d_bump(x,y)
    dbump([x,y],[5,5],[5,5])[1]
end
##
x=range(-10,10, length=10000)
y=copy(x)
gr()
surface(0:1:10,0:1:10,two_d_bump)

##
