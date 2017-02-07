


function plotalpha(data, nbseq, maxposition, m){
	var margin = {top: 30, right: 20, bottom: 70, left: 50},
		width = 800  + margin.left + margin.right,
		height = 300 + margin.top + margin.bottom;
 
	console.log("nbseq:", nbseq);
 
	//Create the Scale we will use for the Axis
	var xscale = d3.scale.linear()
				.domain([0, maxposition])
				.range([0, width]);
        
	var yscale = d3.scale.ordinal()
				.rangeRoundBands([0, height], .0)
				.domain(m);
	
	var xAxisBottom = d3.svg.axis()
				.scale(xscale)
				.orient("bottom")
				.tickFormat(d3.format(".2s"));

	var xAxisTop = d3.svg.axis()
				.scale(xscale)
				.orient("top")
				.tickFormat(d3.format(".2s"));
				
	var zoom = d3.behavior.zoom()
				 .scaleExtent([1, 1000])
				// .x([10])
				 .on("zoom", zoomed);

	//~ var yAxis = d3.svg.axis()
				//~ .scale(yscale)
				//~ .orient("left")
				//~ .tickFormat(d3.format(".2s"));
				
	var alphaplot = d3.select(".alphaplot")
				.attr("width", width + margin.left + margin.right)
				.attr("height", height + margin.top + margin.bottom)
			.append("g")
				.attr("stroke-width",1)				
				.attr("stroke","red")
				.attr("transform", "translate(" + margin.left + "," + margin.top + ")")
				.call(zoom);

	alphaplot.append("g")
//				.attr("stroke-width",1)				
				.attr("class", "x-axis top")
				.attr("transform", "translate(0," + 0 + ")")
//				.attr("stroke","red")
//              .style("font-size","30px")
				.call(xAxisTop);


	alphaplot.append("g")
				.attr("stroke-width",1)				
				.attr("class", "x-axis bottom")
				.attr("transform", "translate(0," + height + ")")
//				.attr("stroke","blue")
//				.attr("fill", "red")
//				.attr("font-size","60px")
				.call(xAxisBottom);

	//~ d3.selectAll("line")
 				//~ .attr("stroke-width",1)				
				//~ .style("fill","red")
                //~ .style("stroke","blue");
                //~ 
                
                
                
                
	//~ d3.selectAll("text")
                //~ .style("font-size","20px")
                //~ .style("stroke", "black")
                //~ .style("fill", "black");


	//~ alphaplot.append("g")
				//~ .attr("class", "y axis")
				//~ .call(yAxis);

	
	//~ circle = alphaplot.selectAll("circle")
				//~ .data(data)
	//~ circle.enter().append("circle")
				//~ .attr("cy", function(d) { return yscale(d.seq); })
				//~ .attr("cx", function(d) { return xscale(d.begin); })
				//~ .attr("r",  function(d) { return xscale(30) });	
	
	arrows = alphaplot.selectAll("polygon")
				.data(data)
			.enter().append("polygon")
				.attr("points",function(d) {
									if (d.strand == 0) {
										vpadding =  ((yscale(d.seq)+yscale(d.seq-1))/2 - yscale(d.seq))/2
										up  =       yscale(d.seq) + vpadding 
										vup =       yscale(d.seq) + 2 * vpadding
										down =      yscale(d.seq) - vpadding 
										vdown =     yscale(d.seq) - 2 * vpadding 
										vmiddle =     yscale(d.seq)
										left =      xscale(d.begin)
										right =     xscale(d.end)
										hpadding =  (xscale(d.end)-xscale(d.begin)) / 3
										hmiddle =    left + 2 * hpadding
										pts= left    + "," + up      + " " + 
									     hmiddle + "," + up      + " " +
									     hmiddle + "," + vup     + " " +
									     right   + "," + vmiddle + " " +
									     hmiddle + "," + vdown     + " " +
									     hmiddle + "," + down      + " " +
									     left + "," + down
									     }
									else {
										vpadding =  ((yscale(d.seq)+yscale(d.seq-1))/2 - yscale(d.seq))/2
										up  =       yscale(d.seq) + vpadding 
										vup =       yscale(d.seq) + 2 * vpadding
										down =      yscale(d.seq) - vpadding 
										vdown =     yscale(d.seq) - 2 * vpadding 
										vmiddle =     yscale(d.seq)
										left =      xscale(d.begin)
										right =     xscale(d.end)
										hpadding =  (xscale(d.end)-xscale(d.begin)) / 3
										hmiddle =    left + 2 * hpadding
										pts= left   + "," + vmiddle + " " +
									     hmiddle + "," + vup    + " " +
									     hmiddle + "," + up     + " " +
									     right   + "," + up     + " " +
									     right   + "," + down   + " " +
									     hmiddle + "," + down   + " " +
									     hmiddle + "," + vdown  
												
									}	
									
									
									return pts
									})
				.attr("stroke","black")
				.attr("stroke-width",1)				
				.attr("fill", function(d) {if (d.strand == 0) {return "blue"} else {return "red"}});
	//~ circle.enter().append("circle")
				//~ .attr("cy", function(d) { return yscale(d.seq); })
				//~ .attr("cx", function(d) { return xscale(d.begin); })
				//~ .attr("r",  function(d) { return xscale(30) });	
			
				
	function zoomed() {
				console.log("zoom:",  d3.event.translate, d3.event.scale);
				alphaplot.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
				}	
}
