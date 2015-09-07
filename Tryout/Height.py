#to calculate the height of the ball at certain time

InitialVelocity=5
AccelerationOfGravity=9.81
TIME=0.9

PositionOfTheBall=InitialVelocity*TIME-0.5*AccelerationOfGravity*TIME**2

print "At t=%g s, the height of the ball is %.2f m" %(TIME, PositionOfTheBall)

#alternative way to express the printout sentence

#failed code below, wait to be checked out
print 'At t={t:AccelerationOfGravity} s, the height of the ball is {PositionOfTheBall:.2f}m.'.format(t=TIME, PositionOfTheBall=PositionOfTheBall)
