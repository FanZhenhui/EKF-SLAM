function [zhat, H] = ObserveModel(XX, zloc, jind)
    H = zeros(2, length(XX));
    %Estimated location relative to the robot's belief
    muJ = XX(jind:jind+1); 
    delta = [muJ(1) - XX(1); muJ(2) - XX(2)];
    q = delta.' * delta;

    angle = AngleLimit(atan2(delta(2), delta(1)) - XX(3));

%     zloc
    zhat = [sqrt(q); angle]
    H(:,1:3) = [-(delta(1)/sqrt(q)) -(delta(2)/sqrt(q)) 0; (delta(2)/q) -(delta(1)/q) -1];
    H(:,jind:jind+1) = [(delta(1)/sqrt(q)) (delta(2)/sqrt(q)); -(delta(2)/q) (delta(1)/q)];
end            

