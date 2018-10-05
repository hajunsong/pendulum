clc; clear all; close all;

num_body = 6;

m_data = load(sprintf('matlab_body%d.txt', num_body));
c_data = load(sprintf('C_body%d.txt', num_body));

index = 2;

ylabel_text = {'Position [rad]', 'Velocity [rad/s]', 'Acceleration [rad/s^2]'};

for i = 1 : num_body
    figure
    set(gcf,'color',[1,1,1])
    for j = 1 : 3
        subplot(3, 1, j)
        plot(c_data(:,1), c_data(:,index), 'LineWidth', 2)
        hold on
        plot(m_data(:,1), m_data(:,index),'--','LineWidth',2)
        grid on
        xlabel('Time [s]')
        ylabel(ylabel_text(j))
        set(gca,'FontSize',13)
        legend('C','MATLAB')
        index = index + 1;
    end
end