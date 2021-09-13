function my_polarhist(thetap, thetaq, bin, arrow, sz, label, ent)

    colors = cbrewer('div', 'BrBG',11);
    if(~ent)
        cp = colors(8,:);
        cq = colors(4,:);
    else
        cp = colors(10,:);
        cq = colors(2,:);
    end
    
    
    hp = polarhistogram(thetap,bin,'FaceColor',cp,'Normalization','count','FaceAlpha',0.8);
    hold on
    hq = polarhistogram(thetaq,bin,'FaceColor',cq,'Normalization','count','FaceAlpha',0.8);
    
    
    ap = 1.2*max(hp.BinCounts);
    aq = 1.2*max(hp.BinCounts);
    % ap=1;aq=1

    h = polarscatter(thetap,ap*ones(size(thetap)), sz,'MarkerFaceColor',cp,'MarkerEdgeColor','none','MarkerFaceAlpha',0.3);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    sump = mean(ap*exp(1i*thetap));
    h = polarscatter(angle(sump),abs(sump), sz,'MarkerFaceColor',cp,'MarkerEdgeColor','#000000');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    h = polarscatter(thetaq,aq*ones(size(thetaq)), sz,'MarkerFaceColor',cq,'MarkerEdgeColor','none','MarkerFaceAlpha',0.3);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    sumq = mean(ap*exp(1i*thetaq));
    h = polarscatter(angle(sumq),abs(sumq), sz,'MarkerFaceColor',cq,'MarkerEdgeColor','#000000');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    % rlim([0,1])
    Ax = gca; 
    Ax.RTickLabel = []; 
    Ax.ThetaTickLabel = [];
    
    text(pi, ap*1.3, label, 'horizon', 'center', 'rotation', 0,'FontSize',8);
    
    if(arrow)
        %%%Data %%%%
        resultant_direction = angle(sump);
        resultant_length = abs(sump);
        %%%%arrow head %%%%
        arrowhead_length    = resultant_length/8; % arrow head length relative to resultant_length
        num_arrowlines = 100;
        arrowhead_angle = deg2rad(30); % degrees
        %%%%arrow tip coordinates %%%%
        t1 = repmat(resultant_direction,1,num_arrowlines);
        r1 = repmat(resultant_length,1,num_arrowlines);
        %%%%arrow base coordinates %%%%
        b = arrowhead_length.*tan(linspace(0,arrowhead_angle,num_arrowlines/2));
        theta = atan(b./(resultant_length-arrowhead_length));
        pre_t2 = [theta, -theta];
        r2 = (resultant_length-arrowhead_length)./cos(pre_t2);
        t2 = t1(1)+pre_t2;
        %%%%plot %%%%
        figure(1)
        polarplot([t1(1) t1(1)],[0 r1(1)-0.9*arrowhead_length],'Color','#000000','linewidth',1.5)
        hold on
        polarplot([t1; t2],[r1; r2],'Color','#000000','linewidth',0.2)
    end
    

    
end
