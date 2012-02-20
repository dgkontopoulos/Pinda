/***************************/
//@Author: Adrian "yEnS" Mato Gondelle &amp;amp;amp; Ivan Guardado Castro
//@website: www.yensdesign.com
//@email: yensamg@gmail.com
//@license: Feel free to use it, but keep this credits please!
/***************************/

$(document).ready(function(){
	$(".menu > li").click(function(e){
		switch(e.target.id){
			case "DNA":
				//change status &amp;amp;amp; style menu
				$("#DNA").addClass("active");
				$("#Protein").removeClass("active");
				//display selected division, hide others
				$("div.DNA").fadeIn();
				$("div.Protein").css("display", "none");
			break;
			case "Protein":
				//change status &amp;amp;amp; style menu
				$("#DNA").removeClass("active");
				$("#Protein").addClass("active");
				//display selected division, hide others
				$("div.Protein").fadeIn();
				$("div.DNA").css("display", "none");
			break;
		}
		//alert(e.target.id);
		return false;
	});
});

