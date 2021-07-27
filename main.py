from flask import Flask, render_template, request, send_file
import nozzleCode

app = Flask(__name__)

@app.route("/", methods=["POST","GET"])
def login():
    if request.method == "POST":
        nozzleCode.generateNozzle(
            float(request.form["chamberPressure"]),
            float(request.form["chamberTemperature"]),
            float(request.form["thrust"]),
            float(request.form["mdot"]),
            float(request.form["altitude"]),
            float(request.form["gamma"]), 
            355,
            float(request.form["throatRadius"]),
            float(request.form["chamberRadius"]))
        return render_template("index.html", nozzleLoaded=1)
    else:
        return render_template("index.html", nozzleLoaded=0)
        
@app.route("/generatedNozzle.stl")
def returnNozzle():
    return send_file('generatedNozzle.stl')

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8080, debug=True)

    