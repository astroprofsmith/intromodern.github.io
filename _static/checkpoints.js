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
            document.getElementById('Answer_q2').textContent = "No, although a hotter star would be brighter, the speed of the stars does not affect their temperature.";
	}
	if (document.getElementById('q2a4').checked) {
            document.getElementById('Answer_q2').textContent = "No, the speed of the star alone would not affect the brightness, but it does affect the arrival time of the photons!";
	}
    }
}
