function Out=GetInterIntraEventResiduals(input)
%This class is used to get observed ground motion residuals and 
%seperate out into inter and intra event terms. There is also the potential
%to consider spatial correlation (but not implemented as a Aug 2012) 
%in the inter and intra event terms

%Input variables:
%input - a structure with the following elements
%       .imObservations - the Im values of the observed ground motions
%       .imMedian - the predicted median ground motion
%       .totalSigma - the predicted total standard deviation
%       .interSigma - the inter-event sigma
%       .intraSigma - the intra-event sigma
%       .eventIds - the event ID numbers

%Output variables:
%Out - a structure with the following elements
%       Out.totalResiduals 
%       Out.totalResidualsNormalized 
%       Out.interEventResiduals 
%       Out.interEventResidualsNormalized 
%       Out.intraEventResiduals 
%       Out.intraEventResidualsNormalized 
%       Out.eventData = eventData - a structure which contains the number
%       of different earthquakes, "numEq", the unique "EqId", the number
%       of Gm per Eq, "eventNumGms", and the inter-event sigma "interSigma" 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make local copies of the structured variables

imObservations = input.imObservations;
eventIds = input.eventIds;
median = input.imMedian;
totalSigma = input.totalSigma;
interSigma = input.interSigma;
intraSigma = input.intraSigma;

numObs=length(imObservations);

%now begin the computations

%get the total residuals
totalResiduals = (log(imObservations) - log(median));
totalResidualsNormalized = totalResiduals./totalSigma;
        
%now get the conditional covariance matrix
condCovMatrix = diag(intraSigma.^2);     

%set the event data 
eventData = setEventData(input);

%now get the inter intra residuals
numEqs = eventData.numEqs;

eta=zeros(numEqs,1);
for i=1:numEqs
    eventNumGms=eventData.eventNumGms(i);
    interEventSigma = eventData.interSigma(i);
    indexmin=sum(eventData.eventNumGms(1:i-1))+1;
    indexmax=sum(eventData.eventNumGms(1:i));

    Ccsub=zeros(indexmax-indexmin+1,indexmax-indexmin+1);
    Ccsub(:,:)=condCovMatrix(indexmin:indexmax,indexmin:indexmax);
    totalResidSub=zeros(indexmax-indexmin+1,1);
    totalResidSub=totalResiduals(indexmin:indexmax,1);
    %get the inter-event residual, eta_i
    oneN = ones(eventNumGms,1);
    eta(i,1) = (oneN'*(Ccsub\totalResidSub))/(1/interEventSigma^2 + oneN'*(Ccsub\oneN));  %note C\ is same as inv(C)
end
%now assign the inter-event term to the eventData structure
eventData.interEventResiduals=eta;
eventData.interEventResidualNormalized=eta./eventData.interSigma;

%now need to assign the inter-event residuals to the individual
%ground motions and get the intra-event residual
for i=1:numObs
    %get the EqId of the ground motion
    gmEqId=eventIds(i);
    %now find the inter-event residual for this earthquake Id
    for j=1:numEqs
        eqId=eventData.eventId(j);
        if gmEqId==eqId
            %now set inter-event residual
            interEventResiduals(i,1)=eventData.interEventResiduals(j);
        end
    end
end
%normalize all the interEvent Residuals
interEventResidualsNormalized=interEventResiduals./interSigma;
%now get the intra-event residual
intraEventResiduals = totalResiduals - interEventResiduals;
%normalize all the intraEvent Residuals
intraEventResidualsNormalized=intraEventResiduals./intraSigma;

%store
Out.totalResiduals = totalResiduals;
Out.totalResidualsNormalized = totalResidualsNormalized;
Out.interEventResiduals = interEventResiduals;
Out.interEventResidualsNormalized = interEventResidualsNormalized;
Out.intraEventResiduals = intraEventResiduals;
Out.intraEventResidualsNormalized = intraEventResidualsNormalized;
Out.eventData = eventData;                  
        
%end of function
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eventData = setEventData(input)
    
    eventIds=input.eventIds; 
    nObs=length(input.imObservations);

    numEq=0;    %initialise
    eventData.eventId=0;
    eventData.eventNumGms=0;

    for i=1:nObs
        foundEqId=0;
        for k=1:numEq
            if eventData.eventId(k)==eventIds(i);
                eventData.eventNumGms(k)=eventData.eventNumGms(k)+1;  %number of events
                foundEqId=1;
                break
            end
        end

        if foundEqId==0
            numEq=numEq+1;
            eventData.eventId(numEq)=eventIds(i);
            eventData.eventNumGms(numEq)=1;  %number of events
        end
    end

    eventData.numEqs=numEq;

    GmId=0;
    for i=1:numEq
        %get the last ground motion which is for this earthquake
        GmId=GmId+eventData.eventNumGms(i);
        interEventSigma=input.interSigma(GmId);
        eventData.interSigma(i,1)=interEventSigma;
    end

    return;
%end of function
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
