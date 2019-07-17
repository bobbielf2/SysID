function U = QoIs(samps,qoi)

type = qoi.type;
switch type
    case 'SizeDistr'
        U = sizeDistribution(samps,qoi);
    case 'MatrixComp'
        U = matrixComposition(samps,qoi);
    case 'ThresholdArea'
        U = thresholdArea(samps,qoi);
    case 'Downsampled'
        U = samps(qoi.xsampind,qoi.tsampind);
end

end