function pale_mycolor = pale_RGB_20180131(mycolor,tint_factor)

pale_mycolor(1) = mycolor(1) + (1 - mycolor(1)) * tint_factor;
pale_mycolor(2) = mycolor(2) + (1 - mycolor(2)) * tint_factor;
pale_mycolor(3) = mycolor(3) + (1 - mycolor(3)) * tint_factor;