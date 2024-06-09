function setup(htmlComponent) 
{
    document.getElementById("btn").addEventListener("click", function(event) 
    {
        htmlComponent.Data = 1234; // dummy value to fire event on matlab side
        htmlComponent.sendEventToMATLAB("ButtonClicked",123);
        
    });
}