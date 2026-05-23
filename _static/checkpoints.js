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
}
