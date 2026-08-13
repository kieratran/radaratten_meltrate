% --- Piecewise linear fit: find best breakpoint in search window ---
% Returns best index (in cropped vector), confidence, and fit coefficients
function [bp, conf, ca, cb] = piecewise_fit(d, p, lo, hi, layer)

mid = length(d)/2; % half of the ice thickness is the threshold
switch layer
    case "firn"
        sse = nan(hi, 1);
        for k = lo:5:hi
            r1  = p(1:k)   - polyval(polyfit(d(1:k),   p(1:k),   1), d(1:k));
            r2  = p(k:mid) - polyval(polyfit(d(k:mid), p(k:mid), 1), d(k:mid));
            sse(k) = sum(r1.^2) + sum(r2.^2);
        end
        [sse_best, bp] = min(sse(lo:hi));
        bp   = bp + lo - 1;
        conf = sse_best / nanmean(sse(lo:hi));
        ca   = polyfit(d(1:bp),   p(1:bp),   1);
        cb   = polyfit(d(bp:end), p(bp:end), 1);
    case "saline"
        sse = nan(hi, 1);
        for k = lo:5:hi
            r1  = p(mid:k)   - polyval(polyfit(d(mid:k),   p(mid:k),   1), d( mid:k));
            r2  = p(k:end) - polyval(polyfit(d(k:end), p(k:end), 1), d(k:end));
            sse(k) = sum(r1.^2) + sum(r2.^2);
        end
        [sse_best, bp] = min(sse(lo:hi));
        bp   = bp + lo - 1;
        conf = sse_best / nanmean(sse(lo:hi));
        ca   = polyfit(d(1:bp),   p(1:bp),   1);
        cb   = polyfit(d(bp:end), p(bp:end), 1);
end
end