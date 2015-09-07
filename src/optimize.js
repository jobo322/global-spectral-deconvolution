/**
 * Created by acastillo on 9/6/15.
 */
var Opt = require("ml-optimize-lorentzian");

function sampleFunction(from, to, x, y, lastIndex){
    //console.log(from+" "+to);
    //console.log(lastIndex+" "+x[lastIndex[0]]);
    var nbPoints = x.length;
    var sampleX = [];
    var sampleY = [];
    var direction = Math.sign(x[1]-x[0]);//Direction of the derivative
    var delta = (to-from)/2;
    var mid = (from+to)/2;
    var stop = false;
    var index = lastIndex[0];
    while(!stop&&index<nbPoints){
        if(Math.abs(x[index]-mid)<=delta){
            sampleX.push(x[index]);
            sampleY.push(y[index]);
            index++;

        }
        //It is outside the range.
        else{

            if(Math.sign(mid-x[index])==direction){
                //We'll reach the mid going in the current direction
                index++;
            }
            else{
                //There is not more peaks in the current range
                stop=true;
            }
        }
        //console.log(sampleX);
    }
    lastIndex[0]=index;
    return [sampleX, sampleY];
}

function optimizePeaks(peakList,x,y,n){
    var i, j, lastIndex=[0];
    var groups = groupPeaks(peakList,n);
    var result = [];
    //console.log(x[0]+" "+x[1]);
    for(i=0;i<groups.length;i++){
        //console.log(peakList[i]);
        var peaks = groups[i].group;
        if(peaks.length>1){
            //Multiple peaks
            console.log("Pending group of overlaped peaks "+peaks.length);
        }
        else{
            //Single peak
            peaks = peaks[0];
            var sampling = sampleFunction(peaks.x-n*peaks.width,
                peaks.x+n*peaks.width,x,y,lastIndex);
            //console.log(sampling);
            if(sampling[0].length>5){
                var error = peaks.width/10000;
                var opts = [  3,    100, error, error, error, error*10, error*10,    11,    9,        1 ];
                //var gauss = Opt.optimizeSingleGaussian(sampling[0], sampling[1], opts, peaks);
                var gauss = Opt.optimizeSingleGaussian2(sampling[0], sampling[1], opts, peaks);
                //console.log(gauss.X2);
                gauss = gauss.p;
                //console.log("Before");
                //console.log(peakList[i]);
                result.push({x:gauss[0],y:gauss[2],width:gauss[1]*2.35482}); // From https://en.wikipedia.org/wiki/Gaussian_function#Properties}
                //console.log("After");
                //console.log(result);
            }
        }

    }
    return result;
}

function groupPeaks(peakList,nL){
    var group = [];
    var groups = [];
    var i, j;
    var limits = [peakList[0].x,nL*peakList[0].width];
    var upperLimit, lowerLimit;
    //Merge forward
    for(i=0;i<peakList.length;i++){
        //If the 2 things overlaps
        if(Math.abs(peakList[i].x-limits[0])<(nL*peakList[i].width+limits[1])){
            //Add the peak to the group
            group.push(peakList[i]);
            //Update the group limits
            upperLimit = limits[0]+limits[1];
            if(peakList[i].x+nL*peakList[i].width>upperLimit){
                upperLimit = peakList[i].x+nL*peakList[i].width;
            }
            lowerLimit = limits[0]-limits[1];
            if(peakList[i].x-nL*peakList[i].width<lowerLimit){
                lowerLimit = peakList[i].x-nL*peakList[i].width;
            }
            limits = [(upperLimit+lowerLimit)/2,Math.abs(upperLimit-lowerLimit)/2];

        }
        else{
            groups.push({limits:limits,group:group});
            //var optmimalPeak = fitSpectrum(group,limits,spectrum);
            group=[peakList[i]];
            limits = [peakList[i].x,nL*peakList[i].width];
        }
    }
    groups.push({limits:limits,group:group});
    //Merge backward
    for(i =groups.length-2;i>=0;i--){
        //The groups overlaps
        if(Math.abs(groups[i].limits[0]-groups[i+1].limits[0])<
            (groups[i].limits[1]+groups[i+1].limits[1])/2){
            for(j=0;j<groups[i+1].group.length;j++){
                groups[i].group.push(groups[i+1].group[j]);
            }
            upperLimit = groups[i].limits[0]+groups[i].limits[1];
            if(groups[i+1].limits[0]+groups[i+1].limits[1]>upperLimit){
                upperLimit = groups[i+1].limits[0]+groups[i+1].limits[1];
            }
            lowerLimit = groups[i].limits[0]-groups[i].limits[1];
            if(groups[i+1].limits[0]-groups[i+1].limits[1]<lowerLimit){
                lowerLimit = groups[i+1].limits[0]-groups[i+1].limits[1];
            }
            //console.log(limits);
            groups[i].limits = [(upperLimit+lowerLimit)/2,Math.abs(upperLimit-lowerLimit)/2];

            groups.splice(i+1,1);
        }
    }
    return groups;
}



module.exports=optimizePeaks;
