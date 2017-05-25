'use strict';

var Opt = require('ml-optimize-lorentzian');

const defaultOptions = {
    widthFactor: 4,
    damping: 0.01,
    errorTolerance: 10e-6,
    maxIterations: 500
};

function optimizePeaks(peakList, x, y, options = {}) {
    options = Object.assign({}, defaultOptions, options);

    var lastIndex = [0];
    var result = [];
    var sampling;

    var factor, optimizeFunctionSum, optimizeFunctionSingle;
    switch (options.functionName) {
        case 'lorentzian':
            factor = 1;
            optimizeFunctionSum = Opt.optimizeLorentzianSum;
            optimizeFunctionSingle = Opt.optimizeSingleLorentzian;
            break;
        default:
            factor = 1.17741; //From https://en.wikipedia.org/wiki/Gaussian_function#Properties
            optimizeFunctionSum = Opt.optimizeGaussianSum;
            optimizeFunctionSingle = Opt.optimizeSingleGaussian;
    }

    var groups = groupPeaks(peakList, options.widthFactor);

    for (let i = 0, l = groups.length; i < l; i++) {
        var peaks = groups[i].group;
        if (peaks.length > 1) {
            sampling = sampleFunction(groups[i].limits[0] - groups[i].limits[1], groups[i].limits[0] + groups[i].limits[1], x, y, lastIndex);
            if (sampling.x.length > 5) {
                var optPeaks = optimizeFunctionSum(sampling, peaks, options);
                for (let j = 0; j < optPeaks.length; j += 3) {
                    result.push({x: optPeaks[j], y: optPeaks[j + 1], width: optPeaks[j + 2] * factor});
                }
            }
        } else {
            peaks = peaks[0];
            sampling = sampleFunction(peaks.x - options.widthFactor * peaks.width,
                peaks.x + options.widthFactor * peaks.width, x, y, lastIndex);
            if (sampling.x.length > 5) {
                var optPeak = optimizeFunctionSum(sampling, [peaks], options);
                result.push({x: optPeak[0], y: optPeak[1], width: optPeak[2] * factor});
            }
        }
    }
    return result;
}

function sampleFunction(from, to, x, y, lastIndex) {
    var nbPoints = x.length;
    var sampleX = [];
    var sampleY = [];
    var direction = Math.sign(x[1] - x[0]);//Direction of the derivative
    if (direction === -1) {
        lastIndex[0] = x.length - 1;
    }
    var delta = Math.abs(to - from) / 2;
    var mid = (from + to) / 2;
    var stop = false;
    var index = lastIndex[0];
    while (!stop && index < nbPoints && index >= 0) {
        if (Math.abs(x[index] - mid) <= delta) {
            sampleX.push(x[index]);
            sampleY.push(y[index]);
            index += direction;
        } else {
            if (Math.sign(mid - x[index]) === 1) {
                index += direction;
            } else {
                stop = true;
            }
        }
    }
    lastIndex[0] = index;
    return {x: sampleX, y: sampleY};
}

function groupPeaks(peakList, nL) {
    var group = [];
    var groups = [];
    var i, j;
    var limits = [peakList[0].x, nL * peakList[0].width];
    var upperLimit, lowerLimit;
    //Merge forward
    for (i = 0; i < peakList.length; i++) {
        //If the 2 things overlaps
        if (Math.abs(peakList[i].x - limits[0]) < (nL * peakList[i].width + limits[1])) {
            //Add the peak to the group
            group.push(peakList[i]);
            //Update the group limits
            upperLimit = limits[0] + limits[1];
            if (peakList[i].x + nL * peakList[i].width > upperLimit) {
                upperLimit = peakList[i].x + nL * peakList[i].width;
            }
            lowerLimit = limits[0] - limits[1];
            if (peakList[i].x - nL * peakList[i].width < lowerLimit) {
                lowerLimit = peakList[i].x - nL * peakList[i].width;
            }
            limits = [(upperLimit + lowerLimit) / 2, Math.abs(upperLimit - lowerLimit) / 2];

        } else {
            groups.push({limits: limits, group: group});
            group = [peakList[i]];
            limits = [peakList[i].x, nL * peakList[i].width];
        }
    }
    groups.push({limits: limits, group: group});
    for (i = groups.length - 2; i >= 0; i--) {
        if (Math.abs(groups[i].limits[0] - groups[i + 1].limits[0]) <
            (groups[i].limits[1] + groups[i + 1].limits[1]) / 2) {
            for (j = 0; j < groups[i + 1].group.length; j++) {
                groups[i].group.push(groups[i + 1].group[j]);
            }
            upperLimit = groups[i].limits[0] + groups[i].limits[1];
            if (groups[i + 1].limits[0] + groups[i + 1].limits[1] > upperLimit) {
                upperLimit = groups[i + 1].limits[0] + groups[i + 1].limits[1];
            }
            lowerLimit = groups[i].limits[0] - groups[i].limits[1];
            if (groups[i + 1].limits[0] - groups[i + 1].limits[1] < lowerLimit) {
                lowerLimit = groups[i + 1].limits[0] - groups[i + 1].limits[1];
            }
            //console.log(limits);
            groups[i].limits = [(upperLimit + lowerLimit) / 2, Math.abs(upperLimit - lowerLimit) / 2];

            groups.splice(i + 1, 1);
        }
    }
    return groups;
}
/**
 * This function try to join the peaks that seems to belong to a broad signal in a single broad peak.
 * @param peakList
 * @param options
 */
function joinBroadPeaks(peakList, options) {
    var width = options.width;
    var broadLines = [];

    var max = 0,
        maxI = 0,
        count = 1;
    for (let i = peakList.length - 1; i >= 0; i--) {
        if (peakList[i].soft) {
            broadLines.push(peakList.splice(i, 1)[0]);
        }
    }
    broadLines.push({x: Number.MAX_VALUE});
    var candidates = {x: [broadLines[0].x], y: [broadLines[0].y]};
    var indexes = [0];
    for (let i = 1; i < broadLines.length; i++) {
        if (Math.abs(broadLines[i - 1].x - broadLines[i].x) < width) {
            candidates.x.push(broadLines[i].x);
            candidates.y.push(broadLines[i].y);
            if (broadLines[i].y > max) {
                max = broadLines[i].y;
                maxI = i;
            }
            indexes.push(i);
            count++;
        } else {
            if (count > 2) {
                var fitted = Opt.optimizeSingleLorentzian(candidates,
                    {x: broadLines[maxI].x, y: max, width: Math.abs(candidates.x[0] - candidates.x[candidates.length - 1])});
                peakList.push({x: fitted[0], y: fitted[1], width: fitted[2], soft: false});

            } else {
                //Put back the candidates to the signals list
                indexes.map(function (index) {
                    peakList.push(broadLines[index]);
                });
            }
            candidates = {x:[broadLines[i].x], y:[broadLines[i].y]};
            indexes = [i];
            max = broadLines[i].y;
            maxI = i;
            count = 1;
        }
    }

    peakList.sort(function (a, b) {
        return a.x - b.x;
    });

    return peakList;

}

module.exports = {optimizePeaks: optimizePeaks, joinBroadPeaks: joinBroadPeaks};

