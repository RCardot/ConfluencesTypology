function idxClosest=ClosestPoint(PointCloud,ReferencePoint)



for i=1:length(ReferencePoint)
    squareDiff=@(a,b) (a-b).^2;
    SquareddistOnComponents(:,i)=bsxfun(squareDiff,PointCloud(:,i),ReferencePoint(i));
end

distToReferencePoint=nan(size(SquareddistOnComponents,1),1);
for j=1:size(SquareddistOnComponents,1)
    distToReferencePoint(j)=sqrt(sum(SquareddistOnComponents(j,:)));
end

% disp(distToReferencePoint)
    idxClosest=distToReferencePoint==min(distToReferencePoint);
end