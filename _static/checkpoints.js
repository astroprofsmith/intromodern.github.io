function validate_ans(myq) {
    // What is the paramter being changed in the MM experiment
    if (myq == 'q1') {
	if (document.getElementById('q1a1').checked) {
            document.getElementById('Answer_q1').textContent = "You might be confusing this experiment with other interferometery experiments, where one varies the number of wavelengths along each path by changing the length of the path.  In this experiment, it's important that the paths stay the same length.";
	}
	if (document.getElementById('q1a2').checked) {
            document.getElementById('Answer_q1').textContent = "If there were an ether, the wavelength of the light would change in a direction-dependent way, but that would not be the variable that is being changed by the scientists running the experiment.";
	}
	if (document.getElementById('q1a3').checked) {
            document.getElementById('Answer_q1').textContent = "Yes!  By rotating the apparatus, one ensures that if there were an ether, the light would sometimes be going with the ether and sometimes against it.";
	}
	if (document.getElementById('q1a4').checked) {
            document.getElementById('Answer_q1').textContent = "Technically, one cannot change the speed of the ether.  If there were an ether, and one did this experiment at different times, the motion of the Earth through the ether would be different, so one could say that one has changed the relative speed of the apparatus and the ether, but that is not technically what the scientists are changing.";
	}
    }




    if (myq == 'q2') {
	// Why is the light curve brighter?
	if (document.getElementById('q2a1').checked) {
            document.getElementById('Answer_q2').textContent = "Yes!  If light got a speed boost from the motion of the stars, by the time the light got here, we would see changes in brightness depending on how the stars were moving when the light was emitted.  We don't see that.";
	}
	if (document.getElementById('q2a2').checked) {
            document.getElementById('Answer_q2').textContent = "No -- that does explain why one of the eclipse dips is deeper than the other, but not why the brightness of BOTH stars together would get brighter.";
	}
	if (document.getElementById('q2a3').checked) {
            document.getElementById('Answer_q2').textContent = "No, although a hotter star would be brighter, the speed of the stars does not affect their temperature.  The speeds involved are nowhere near what you would need for relativitisic Doppler boosting to make approaching star seem brighter and hotter than the receding.";
	}
	if (document.getElementById('q2a4').checked) {
            document.getElementById('Answer_q2').textContent = "No, the speed of the star alone would not affect the brightness, but it does affect the arrival time of the photons!";
	}
    }

    if (myq == 'q3') {
	// Clock reading or duration?
	if (document.getElementById('q3a1').checked) {
            document.getElementById('Answer_q3').textContent = "The 08:32 would nominally be a clock reading, although it could of course be the duration since the most recent noon or midnight.";
	}
	if (document.getElementById('q3a2').checked) {
            document.getElementById('Answer_q3').textContent = "This would be the difference in two clock readings (stop and start) so therefore a duration.  If the start time were zero, it could also be a clock reading.";
	}
	if (document.getElementById('q3a3').checked) {
            document.getElementById('Answer_q3').textContent = "This most likely a duration -- a second event happened one day after a first event.  However, it is of course possible that the first clock reading was zero, and therefore the clock reading at the end of the duration would also be one day.";
	}
	if (document.getElementById('q3a4').checked) {
            document.getElementById('Answer_q3').textContent = "This would be a clock reading for a sundial.  It could also represent a duration, if you measured how much the angle changed since some earlier time in the day.";
	}
    }

    if (myq == 'q4') {
	// Ruler reading, distance, or displacement?
	if (document.getElementById('q4a1').checked) {
            document.getElementById('Answer_q4').textContent = "This is a little tricky, because of the word 'over'.  If there were a word like 'right' or 'left' included, it would be a displacement. Since no direction is specified, it would have to be a distance of two seats, with the direction left ambiguous.";
	}
	if (document.getElementById('q4a2').checked) {
            document.getElementById('Answer_q4').textContent = "This would be a ruler reading -- a mark on a specific place on the ruler.  It could also be a distance of 1.98 m from the end of the tape.";
	}
	if (document.getElementById('q4a3').checked) {
            document.getElementById('Answer_q4').textContent = "This is a distance.  The difference in ruler readings between the top and the bottom of the dog.";
	}
	if (document.getElementById('q4a4').checked) {
            document.getElementById('Answer_q4').textContent = "This is a displacement: distance and direction.";
	}
    }


    if (myq == 'q5') {
	// Inertial reference frames
	if (document.getElementById('q5a1').checked) {
            document.getElementById('Answer_q5').textContent = "If your room is in a standard building on stable ground, it can be considered an inertial reference frame.";
	}
	if (document.getElementById('q5a2').checked) {
            document.getElementById('Answer_q5').textContent = "If the elevator is in free-fall, objects will drift with apparaent constant velocities in the elevator.  It does act like an inertial reference frame.";
	}
	if (document.getElementById('q5a3').checked) {
            document.getElementById('Answer_q5').textContent = "As long as the spaceship doesn't change its velocity, it will act as an inertial reference frame.  If it rotates, however, an apparent force inside the ship will seem to push objects toward the hull, making an artificial gravity.  That would not be an inertial frame.";
	}
	if (document.getElementById('q5a4').checked) {
            document.getElementById('Answer_q5').textContent = "From the inside of the plane, it would seem like mysterious forces were pushing the people and objects around the cabin.  This is not an inertial reference frame.";
	}
    }

    if (myq == 'q6') {
	// Operational definitions
	if (document.getElementById('q6a1').checked) {
            document.getElementById('Answer_q6').textContent = "No, this is at best a poetic expression, not a definition, and certainly not an operational definition.  There's no procedure or measurement.";
	}
	if (document.getElementById('q6a2').checked) {
            document.getElementById('Answer_q6').textContent = "No, this is more a metaphor, and a confusing one at that.";
	}
	if (document.getElementById('q6a3').checked) {
            document.getElementById('Answer_q6').textContent = "This is operational.  Stop an object from falling and measure the force exerted to do so.  Call that force the weight.";
	}
	if (document.getElementById('q6a4').checked) {
            document.getElementById('Answer_q6').textContent = "This is operational.  Measure the volume of the object, exert the standard pressure, measure the new volume, calculate the percentage change, and call that softness.";
	}
	if (document.getElementById('q6a5').checked) {
            document.getElementById('Answer_q6').textContent = "This is not operational.  There is no procedure and nothing being measured or calculated.";
	}
    }


    if (myq =='q7') {
	// Postulates and implications
	if (document.getElementById('q7a1').checked) {
            document.getElementById('Answer_q7').textContent = "I cannot just know how far away a distant event is.  I would need to either assemble clocks and rulers, or send a signal to you at a known speed, to figure out a remote distance.";
	}
	if (document.getElementById('q7a2').checked) {
            document.getElementById('Answer_q7').textContent = "The catching is an effect, and the throwing is a cause.  Therefore, it is not possible to observe these particular events in the opposite order. (But as you will see, events that are not cause and effect can be in the other order from a different reference frame.)";
	}
	if (document.getElementById('q7a3').checked) {
            document.getElementById('Answer_q7').textContent = "No.  We will see that the runner (you) might measure a different duration for the race on their own stopwatch than the person on the field (me), but you and I need to agree what my stopwatch actually says.";
	}
	if (document.getElementById('q7a4').checked) {
            document.getElementById('Answer_q7').textContent = "Events must still happen in any reference frame.  The time and space coordinates of events may well be different in different frames, but the events have to still happen.";
	}
    }


    if (myq =='q8') {
	//  Half the sped of light
	if (document.getElementById('q8a1').checked) {
            document.getElementById('Answer_q8').textContent = "This is a ludicrously small speed compared to the speed of light.  Much less than half.";
	}
	if (document.getElementById('q8a2').checked) {
            document.getElementById('Answer_q8').textContent = "Yes!  This is half of the speed of light.";
	}
	if (document.getElementById('q8a3').checked) {
            document.getElementById('Answer_q8').textContent = "This is the escape velocity of the Earth.  Still small compared to c.";
	}
	if (document.getElementById('q8a4').checked) {
            document.getElementById('Answer_q8').textContent = "This is the actual speed of light.  beta = 1, not one-half.";
	}
    }

    if (myq =='q9') {
	//  Time dilation
	if (document.getElementById('q9a1').checked) {
            document.getElementById('Answer_q9').textContent = "Adam is at rest with respect to the two events, and will therefore measure the shortest time interval (the proper time interval).  Of the four reference frames, this clock is running slowest.";
	}
	if (document.getElementById('q9a2').checked) {
            document.getElementById('Answer_q9').textContent = "The relative speed here is not zero, so Braden will measure a longer time interval than Adam, but not by much (50 mph is *very* slow, compared to c).";
	}
	if (document.getElementById('q9a3').checked) {
            document.getElementById('Answer_q9').textContent = "Even at 500 km/hr, this is still very, very slow compared to c.  Cathy will measure a (slightly) longer time interval than Adam or Braden, but not the longest.";
	}
	if (document.getElementById('q9a4').checked) {
            document.getElementById('Answer_q9').textContent = "At this speed, Zaphod will measure an interval roughly twice as long as Adam does.  This is the longest duration of the four, and therefore was can say Zaphod's clock is running the fastest.  This is why 'moving clocks run slow' is so confusing. The phrase itself does not make it clear which clock is moving.";
	}
    }


    if (myq =='q10') {
	//  How to set up inverse transform
	if (document.getElementById('q10a1').checked) {
            document.getElementById('Answer_q10').textContent = "No, you don't want x to point the other way.  The relative motion is in the other direction, not which way x is increasing.";
	}
	if (document.getElementById('q10a2').checked) {
            document.getElementById('Answer_q10').textContent = "Yes, this.  You could also switch which of the two frames you call prime -- the starting frame is always on top.";
	}
	if (document.getElementById('q10a3').checked) {
            document.getElementById('Answer_q10').textContent = "No, time axes always point up.";
	}
	if (document.getElementById('q10a4').checked) {
            document.getElementById('Answer_q10').textContent = "No, we keep the same orientation of the axes -- we just want to switch the relative direction of the motion.";
	}
    }

    if (myq =='q11') {
	//  Properties of four vectors
	if (document.getElementById('q11a1').checked) {
            document.getElementById('Answer_q11').textContent = "This is accurate.  All four-vectors will have a negative squared time component.  Even if the time component itself is negative, that minus will become plus in the square, while the i (or the covariant term) will ensure a minus.";
	}
	if (document.getElementById('q11a2').checked) {
            document.getElementById('Answer_q11').textContent = "This is accurate.  This follows from the requirement that c be the same in all reference frames.";
	}
	if (document.getElementById('q11a3').checked) {
            document.getElementById('Answer_q11').textContent = "This is not accurate.  You CAN put an i on the time component, but you do not HAVE to.  It's not really an imaginary number, either -- it's a way of keeping track of the minus sign.";
	}
	if (document.getElementById('q11a4').checked) {
            document.getElementById('Answer_q11').textContent = "This is accurate.  We do not always write down all four components, but they are there.";
	}
    }

    if (myq =='q12') {
	//  Values of the Lorentz Transformation
	if (document.getElementById('q12a1').checked) {
            document.getElementById('Answer_q12').textContent = "No.";
	}
	if (document.getElementById('q12a2').checked) {
            document.getElementById('Answer_q12').textContent = "No.";
	}
	if (document.getElementById('q12a3').checked) {
            document.getElementById('Answer_q12').textContent = "No.";
	}
	if (document.getElementById('q12a4').checked) {
            document.getElementById('Answer_q12').textContent = "Yes!  That's it!";
	}
    }    

    if (myq =='q13') {
	//  What changes and what doesn't under a LT
	if (document.getElementById('q13a1').checked) {
            document.getElementById('Answer_q13').textContent = "Yes, the x-axis value will change.";
	}
	if (document.getElementById('q13a2').checked) {
            document.getElementById('Answer_q13').textContent = "No, the y-axis value will not change.";
	}
	if (document.getElementById('q13a3').checked) {
            document.getElementById('Answer_q13').textContent = "Yes, the speed of an object will change, unless it is already going at the speed of light (c).";
	}
	if (document.getElementById('q13a4').checked) {
            document.getElementById('Answer_q13').textContent = "Yes, the angle will change, unless it is zero or pi.";
	}
	if (document.getElementById('q13a5').checked) {
            document.getElementById('Answer_q13').textContent = "Yes, the displacement will be different, even if the object is stationary -- it won't be stationary in the second frame (plus time dilation will change the time-displacement).";
	}
    }


    if (myq =='q14') {
	//  this is a template to make new questions
	if (document.getElementById('q14a1').checked) {
            document.getElementById('Answer_q14').textContent = "No, in spacetime diagrams, the vertical axis is time, not height.";
	}
	if (document.getElementById('q14a2').checked) {
            document.getElementById('Answer_q14').textContent = "No, higher up means later.";
	}
	if (document.getElementById('q14a3').checked) {
            document.getElementById('Answer_q14').textContent = "Yes! Lower down means earlier.";
	}
	if (document.getElementById('q14a4').checked) {
            document.getElementById('Answer_q14').textContent = "No, the events would have to be at the same height to be at the same time in this diagram.";
	}
    }    
    
    if (myq =='q15') {
	//  this is a template to make new questions
	if (document.getElementById('q15a1').checked) {
            document.getElementById('Answer_q15').textContent = "True.  For light, the sides of the triangle are equal, so the slope is one.";
	}
	if (document.getElementById('q15a2').checked) {
            document.getElementById('Answer_q15').textContent = "False.  EVERY event has its own set of light cones.";
	}
	if (document.getElementById('q15a3').checked) {
            document.getElementById('Answer_q15').textContent = "False.  The DEFINITION comes from calculating with t1 and t2.  We interpret those numbers as lengths of the sides of a triangle in Aristotalean space and time, but those triangle sides are not the definition.";
	}
	if (document.getElementById('q15a4').checked) {
            document.getElementById('Answer_q15').textContent = "True.  Although you could use the convention that the time component (squared) is positive, in the convention I am using in this book, space terms are positive, so if the space displacement is bigger than the time displacement, the difference of the squares will be positive.";
	}
    }    

    if (myq =='q16') {
	//  How to line up the axes
	if (document.getElementById('q16a1').checked) {
            document.getElementById('Answer_q16').textContent = "Neither.  Cyan is to the right, red happens after.";
	}
	if (document.getElementById('q16a2').checked) {
            document.getElementById('Answer_q16').textContent = "Cyan is to the right, red is at the same time.";
	}
	if (document.getElementById('q16a3').checked) {
            document.getElementById('Answer_q16').textContent = "Cyan is in the same place, red happens earlier.";
	}
	if (document.getElementById('q16a4').checked) {
            document.getElementById('Answer_q16').textContent = "Neither.  Cyan is to the left, and red happens earlier.";
	}
    }
    
    if (myq =='q17') {
	//  rapidity and angle
	if (document.getElementById('q17a1').checked) {
            document.getElementById('Answer_q17').textContent = "No, that will not move the axis at all.";
	}
	if (document.getElementById('q17a2').checked) {
            document.getElementById('Answer_q17').textContent = "Yes!  For beta=0.578, the angle will be 30 deg.";
	}
	if (document.getElementById('q17a3').checked) {
            document.getElementById('Answer_q17').textContent = "No, that corresponds to beta=0.856 and a 40 deg angle.";
	}
	if (document.getElementById('q17a4').checked) {
            document.getElementById('Answer_q17').textContent = "No, that would be very close to a 45 deg angle.";
	}
    }

   if (myq =='q18') {
	//  What does not change in a Lorentz transformation?
	if (document.getElementById('q18a1').checked) {
            document.getElementById('Answer_q18').textContent = "This will not change.  It's important to recognize that an observer in a moving frame will measure a different time interval, but must still agree that the observer in the rest frame measures the proper time interval.";
	}
	if (document.getElementById('q18a2').checked) {
            document.getElementById('Answer_q18').textContent = "This will change.";
	}
	if (document.getElementById('q18a3').checked) {
            document.getElementById('Answer_q18').textContent = "Since the displacement in x changes, the total displacement through space must also change.  The displacement through spacetime remains the same (the size of the four-displacement), but the displacements through space and time separately both change.";
	}
	if (document.getElementById('q18a4').checked) {
            document.getElementById('Answer_q18').textContent = "This will change.  This is time dilation.";
	}
    }    


   
    if (myq =='q19') {
	//  True statement about muons
	if (document.getElementById('q19a1').checked) {
            document.getElementById('Answer_q19').textContent = "Yes, this is the distance traveled as measured in the Earth's rest frame.";
	}
	if (document.getElementById('q19a2').checked) {
            document.getElementById('Answer_q19').textContent = "Yes, this is how fast the muons are traveling, and this is how fast the muons would think the Earth is traveling.";
	}
	if (document.getElementById('q19a3').checked) {
            document.getElementById('Answer_q19').textContent = "Yes, this observation has to be valid in all reference frames.";
	}
	if (document.getElementById('q19a4').checked) {
            document.getElementById('Answer_q19').textContent = "No, the lifetime of a particle is always reported in its own rest frame, since it will be longer in any other frame (and in particular, longer in the Earth's rest frame).";
	}
    }    


    if (myq =='q20') {
	//  Why not use dt0?
	if (document.getElementById('q20a1').checked) {
            document.getElementById('Answer_q20').textContent = "This isn't wrong, but is it the only valid answer?";
	}
	if (document.getElementById('q20a2').checked) {
            document.getElementById('Answer_q20').textContent = "Not wrong, but not the only one.";
	}
	if (document.getElementById('q20a3').checked) {
            document.getElementById('Answer_q20').textContent = "Sure, could be, but is there nothing else?";
	}
	if (document.getElementById('q20a4').checked) {
            document.getElementById('Answer_q20').textContent = "Yup!  It's all three of them.";
	}
    }    

 
    if (myq =='q21') {
	//  this is a template to make new questions
	if (document.getElementById('q21a1').checked) {
            document.getElementById('Answer_q21').textContent = "No, the person on the planet would still have to see the alien going faster than you.";
	}
	if (document.getElementById('q21a2').checked) {
            document.getElementById('Answer_q21').textContent = "No, the person on the planet would have to measure a faster speed for the alien than you do.";
	}
	if (document.getElementById('q21a3').checked) {
            document.getElementById('Answer_q21').textContent = "Yes, this is how your speed and the alien's speed would add.";
	}
	if (document.getElementById('q21a4').checked) {
            document.getElementById('Answer_q21').textContent = "No, no physical object with mass, even an alien, can move at c.";
	}
    }    


    if (myq =='q22') {
	//  this is a template to make new questions
	if (document.getElementById('q22a1').checked) {
            document.getElementById('Answer_q22').textContent = "This statement is not wrong.  If an object needs to actually move through the displacement between two events, then the four displacement needs to be timelike, yes, but one can define a four-displacement between any two events.  However, this has nothing to do with the perpendicular direction issue.";
	}
	if (document.getElementById('q22a2').checked) {
            document.getElementById('Answer_q22').textContent = "Yes, all components of the four velocity are affected by time dilation, and this provides a channel by which motion in x can affect motion in y.  This becomes even more extreme when we consider four-acceleration!";
	}
	if (document.getElementById('q22a3').checked) {
            document.getElementById('Answer_q22').textContent = "This is just not true.";
	}
	if (document.getElementById('q22a4').checked) {
            document.getElementById('Answer_q22').textContent = "No, there's a gamma in the time component, which is not a constant.";
	}
    }

    if (myq =='q23') {
	//  Why use gamma in momentum?
	if (document.getElementById('q23a1').checked) {
            document.getElementById('Answer_q23').textContent = "This isn't wrong, but is it the only valid answer?";
	}
	if (document.getElementById('q23a2').checked) {
            document.getElementById('Answer_q23').textContent = "Not wrong, but not the only one.";
	}
	if (document.getElementById('q23a3').checked) {
            document.getElementById('Answer_q23').textContent = "Sure, that makes sense, but is there nothing else?";
	}
	if (document.getElementById('q23a4').checked) {
            document.getElementById('Answer_q23').textContent = "Yup!  It's all three of them.";
	}
    }    

    if (myq =='q24') {
	//  Sun destroys mass into light energy
	if (document.getElementById('q24a1').checked) {
            document.getElementById('Answer_q24').textContent = "Not nearly enough.";
	}
	if (document.getElementById('q24a2').checked) {
            document.getElementById('Answer_q24').textContent = "Yes! To produce that much energy, you would need to destroy this much mass.";
	}
	if (document.getElementById('q24a3').checked) {
            document.getElementById('Answer_q24').textContent = "It would be five billion kg, not tons.";
	}
	if (document.getElementById('q24a4').checked) {
            document.getElementById('Answer_q24').textContent = "Mass is not conserved!";
	}
    }

    if (myq =='q25') {
	//  Acceleration from force
	if (document.getElementById('q25a1').checked) {
            document.getElementById('Answer_q25').textContent = "No, too slow.  You might be ignoring units.";
	}
	if (document.getElementById('q25a2').checked) {
            document.getElementById('Answer_q25').textContent = "No, too slow.";
	}
	if (document.getElementById('q25a3').checked) {
            document.getElementById('Answer_q25').textContent = "Yes, five divided by 0.5 is ten.";
	}
	if (document.getElementById('q25a4').checked) {
            document.getElementById('Answer_q25').textContent = "No, too fast.";
	}
    }        
    
    if (myq =='q26') {
	//  work to KE theorem
	if (document.getElementById('q26a1').checked) {
            document.getElementById('Answer_q26').textContent = "Yes, that integral ignores the presence of the gamma.";
	}
	if (document.getElementById('q26a2').checked) {
            document.getElementById('Answer_q26').textContent = "Yes, you may integrate a vector.";
	}
	if (document.getElementById('q26a3').checked) {
            document.getElementById('Answer_q26').textContent = "The dot product is not an issue here.";
	}
	if (document.getElementById('q26a4').checked) {
            document.getElementById('Answer_q26').textContent = "As long as speeds are small, you can get away with this, but be careful!  This derivation does not have any speed limit!";
	}
    }

    if (myq =='q27') {
	//  Why F or F'?
	if (document.getElementById('q27a1').checked) {
            document.getElementById('Answer_q27').textContent = "No, the location of the prime does change the meaning.";
	}
	if (document.getElementById('q27a2').checked) {
            document.getElementById('Answer_q27').textContent = "No, then it would be upside down.";
	}
	if (document.getElementById('q27a3').checked) {
            document.getElementById('Answer_q27').textContent = "Always possible, but I am pretty sure I got this one right.";
	}
	if (document.getElementById('q27a4').checked) {
            document.getElementById('Answer_q27').textContent = "Yes, it's there on both sides.  It's just that for F, beta is zero, so the fraction is one.";
	}
    }    

    if (myq =='q28') {
	//  this is a template to make new questions
	if (document.getElementById('q28a1').checked) {
            document.getElementById('Answer_q28').textContent = "This statement is correct.";
	}
	if (document.getElementById('q28a2').checked) {
            document.getElementById('Answer_q28').textContent = "This statement is correct.";
	}
	if (document.getElementById('q28a3').checked) {
            document.getElementById('Answer_q28').textContent = "This statement is too broad.  Electric repulsion can overwhelm magnetic attraction, and currents in opposite directions repel each other.";
	}
	if (document.getElementById('q28a4').checked) {
            document.getElementById('Answer_q28').textContent = "This statement is correct.";
	}
    }    

    
    if (myq =='q2X') {
	//  this is a template to make new questions
	if (document.getElementById('q2Xa1').checked) {
            document.getElementById('Answer_q2X').textContent = "";
	}
	if (document.getElementById('q2Xa2').checked) {
            document.getElementById('Answer_q2X').textContent = "";
	}
	if (document.getElementById('q2Xa3').checked) {
            document.getElementById('Answer_q2X').textContent = "";
	}
	if (document.getElementById('q2Xa4').checked) {
            document.getElementById('Answer_q2X').textContent = "";
	}
    }    
}
